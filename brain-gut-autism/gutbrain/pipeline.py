"""End-to-end gut–brain analysis: statistics, integration, prediction, report."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import microbiome as mb
from . import plots
from .data import Cohort
from .integration import sparse_cca
from .models import regress_severity_cv, run_fusion_cv
from .stats import correlate_blocks, group_differences, mediation

SEVERITY_CANDIDATES = ("SRS_total", "ADOS_CSS", "CARS_total", "ABC_total")


@dataclass
class PipelineConfig:
    covariates: list[str] = field(default_factory=lambda: ["age", "sex"])
    mediation: list[tuple[str, str, str]] = field(default_factory=list)
    severity: str | None = None
    min_prevalence: float = 0.1
    n_perm: int = 500
    cv_repeats: int = 3
    seed: int = 0


def _fmt(df: pd.DataFrame, floatfmt: str = "{:.3g}") -> str:
    """Markdown table without the optional `tabulate` dependency."""
    df = df.copy()
    for c in df.columns:
        if pd.api.types.is_float_dtype(df[c]):
            df[c] = df[c].map(lambda v: "" if pd.isna(v) else floatfmt.format(v))
    head = "| " + " | ".join(map(str, df.columns)) + " |"
    sep = "|" + "---|" * len(df.columns)
    rows = ["| " + " | ".join(map(str, r)) + " |" for r in df.itertuples(index=False)]
    return "\n".join([head, sep, *rows])


def _lookup(name: str, frames: dict[str, pd.DataFrame]) -> pd.Series:
    for df in frames.values():
        if df is not None and name in df.columns:
            return df[name]
    raise KeyError(f"'{name}' not found in any table")


def run_pipeline(cohort: Cohort, out_dir: str | Path, cfg: PipelineConfig) -> dict:
    out = Path(out_dir)
    fig_dir = out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    y = cohort.y
    groups = cohort.subjects["group"]
    cov = cohort.covariates([c for c in cfg.covariates if c in cohort.subjects.columns])
    summary: dict = {"cohort": cohort.summary(), "covariates": list(cov.columns)}
    md = [
        "# Gut–brain axis analysis report",
        "",
        "```",
        cohort.summary(),
        "```",
        "",
        f"Covariates regressed out / adjusted for: `{', '.join(cov.columns) or 'none'}`",
        "",
    ]
    md += [f"> note: {n}" for n in cohort.notes]

    # ---- 1. microbiome ---------------------------------------------------------
    counts = mb.filter_taxa(cohort.microbiome, cfg.min_prevalence)
    clr = mb.clr(counts)
    alpha = mb.alpha_diversity(cohort.microbiome)
    alpha_cov = cov.assign(log_depth=np.log(alpha["depth"]))
    alpha_res = group_differences(alpha[["shannon", "observed", "inv_simpson"]], y, alpha_cov)
    plots.alpha_diversity(alpha, groups, fig_dir / "alpha_diversity.png")

    clr_res = mb.residualize(clr, cov)
    dist = mb.aitchison_distance(clr_res)
    perm = mb.permanova(dist, groups.to_numpy(), n_perm=999, seed=cfg.seed)
    coords, var = mb.pcoa(dist)
    plots.ordination(coords, var, groups, "Aitchison PCoA (covariate-adjusted)", fig_dir / "pcoa.png")

    da = group_differences(clr, y, cov)
    da.to_csv(out / "microbiome_differential_abundance.csv", index=False)
    plots.volcano(da, "Differential abundance (CLR)", fig_dir / "volcano_microbiome.png")
    summary["permanova"] = perm
    summary["n_taxa_tested"] = int(clr.shape[1])
    summary["n_taxa_q05"] = int((da["q_value"] < 0.05).sum())

    md += [
        "## 1. Microbiome",
        "",
        f"{counts.shape[1]} of {cohort.microbiome.shape[1]} taxa kept "
        f"(prevalence ≥ {cfg.min_prevalence:.0%}); CLR-transformed.",
        "",
        "**Alpha diversity** (adjusted for covariates and log sequencing depth):",
        "",
        _fmt(alpha_res[["feature", "effect_d", "p_value", "q_value"]]),
        "",
        f"**Beta diversity**: PERMANOVA on covariate-residualized Aitchison distance — "
        f"pseudo-F = {perm['pseudo_F']:.2f}, R² = {perm['R2']:.3f}, p = {perm['p_value']:.3g}",
        "",
        "![alpha](figures/alpha_diversity.png) ![pcoa](figures/pcoa.png)",
        "",
        f"**Differential abundance** (CLR ~ group + covariates, BH-FDR): "
        f"{summary['n_taxa_q05']} taxa with q < 0.05. Top 10:",
        "",
        _fmt(da.head(10)[["feature", "effect_d", "mean_ASD", "mean_TD", "q_value"]]),
        "",
        "![volcano](figures/volcano_microbiome.png)",
        "",
    ]

    blocks = {"microbiome": clr}
    frames = {"subjects": cohort.subjects, "microbiome": clr}

    # ---- 2. brain & metabolites: univariate ------------------------------------
    metab = cohort.metabolites
    if metab is not None and (metab > 0).all().all() and metab.max().max() > 100:
        metab = np.log2(metab)  # raw intensities -> log2
    for name, df in (("brain", cohort.brain), ("metabolites", metab)):
        if df is None:
            continue
        blocks[name] = df
        frames[name] = df
        res = group_differences(df, y, cov)
        res.to_csv(out / f"{name}_group_differences.csv", index=False)
        summary[f"n_{name}_q05"] = int((res["q_value"] < 0.05).sum())
        md += [
            f"## 2{'a' if name == 'brain' else 'b'}. {name.capitalize()} — ASD vs TD",
            "",
            f"{df.shape[0]} subjects × {df.shape[1]} features; "
            f"{summary[f'n_{name}_q05']} with q < 0.05. Top 10:",
            "",
            _fmt(res.head(10)[["feature", "effect_d", "mean_ASD", "mean_TD", "q_value"]]),
            "",
        ]

    # ---- 3. gut <-> brain integration -------------------------------------------
    if cohort.brain is not None:
        both = clr.index.intersection(cohort.brain.index)
        gx = mb.residualize(clr.loc[both], cov)
        bx = mb.residualize(cohort.brain.loc[both], cov)
        scca = sparse_cca(gx, bx, n_perm=cfg.n_perm, seed=cfg.seed)
        plots.scca_mode(scca, groups, fig_dir / "scca_mode1.png")
        scca.x_weights.to_csv(out / "scca_gut_weights.csv")
        scca.y_weights.to_csv(out / "scca_brain_weights.csv")
        summary["scca"] = {
            "n": len(both),
            "in_sample_r": scca.in_sample_r,
            "cv_r": scca.cv_r,
            "perm_p": scca.perm_p,
        }
        top = lambda w: ", ".join(  # noqa: E731
            f"{k} ({v:+.2f})"
            for k, v in w["mode1"][w["mode1"] != 0].sort_values(key=abs, ascending=False).head(6).items()
        )
        md += [
            "## 3. Gut ↔ brain integration (sparse CCA)",
            "",
            f"n = {len(both)} children with both blocks; both residualized on covariates. "
            f"Mode 1: in-sample r = {scca.in_sample_r[0]:.2f}, cross-validated r = "
            f"{scca.cv_r:.2f}, permutation p = {scca.perm_p:.3g} ({cfg.n_perm} permutations).",
            "",
            f"- Gut side: {top(scca.x_weights)}",
            f"- Brain side: {top(scca.y_weights)}",
            "",
            "![scca](figures/scca_mode1.png)",
            "",
        ]
    if metab is not None:
        both = clr.index.intersection(metab.index)
        da_top = da.head(15)["feature"].tolist()
        corr = correlate_blocks(
            mb.residualize(clr.loc[both, da_top], cov), mb.residualize(metab.loc[both], cov)
        )
        corr.to_csv(out / "taxa_metabolite_correlations.csv", index=False)
        md += [
            "### Taxa ↔ metabolite links (top differential taxa, Spearman, residualized)",
            "",
            _fmt(corr.head(10)),
            "",
        ]

    # ---- 4. mediation ------------------------------------------------------------
    if cfg.mediation:
        md += [
            "## 4. Mediation (hypothesis-driven)",
            "",
            "Specified a priori; adjusted for covariates and diagnosis. "
            "Cross-sectional mediation is compatible with, not proof of, a causal chain.",
            "",
        ]
        rows = []
        med_cov = cov.assign(ASD=y)
        for i, (xn, mn, yn) in enumerate(cfg.mediation):
            res = mediation(
                _lookup(xn, frames), _lookup(mn, frames), _lookup(yn, frames),
                med_cov, seed=cfg.seed,
            )
            rows.append({"X": xn, "M": mn, "Y": yn, **res})
            plots.mediation_diagram(res, xn, mn, yn, fig_dir / f"mediation_{i + 1}.png")
            md.append(f"![mediation {i + 1}](figures/mediation_{i + 1}.png)")
        med_df = pd.DataFrame(rows)
        med_df.to_csv(out / "mediation.csv", index=False)
        summary["mediation"] = rows
        md += [
            "",
            _fmt(med_df[["X", "M", "Y", "n", "a_path", "b_path", "indirect_ab",
                         "indirect_ci_low", "indirect_ci_high", "indirect_p_boot"]]),
            "",
        ]

    # ---- 5. prediction ----------------------------------------------------------
    fusion = run_fusion_cv(blocks, y, cov, n_repeats=cfg.cv_repeats, seed=cfg.seed)
    auc = fusion.summary()
    auc.to_csv(out / "classification_auc.csv")
    fusion.oof_proba.to_csv(out / "classification_oof_probabilities.csv")
    plots.model_auc(auc, fig_dir / "model_auc.png")
    summary["classification"] = auc.round(4).to_dict(orient="index")
    md += [
        "## 5. Prediction of diagnosis (nested CV)",
        "",
        f"Outer 5-fold × {cfg.cv_repeats} repeats, inner 5-fold tuning; covariates "
        "regressed out inside each training fold. Late fusion stacks the per-block "
        "logits, so children with a missing block are still scored. "
        "`AUC_complete_cases` scores every model on the same children (those with "
        "all blocks) for a like-for-like comparison.",
        "",
        _fmt(auc.reset_index(names="model")),
        "",
        "Mean late-fusion weight per block: "
        + ", ".join(f"{k} = {v:.2f}" for k, v in fusion.meta_weights.mean().items()),
        "",
        "![auc](figures/model_auc.png)",
        "",
    ]
    for m, coef in fusion.coefficients.items():
        coef.to_csv(out / f"classifier_coefficients_{m}.csv")
        md += [
            f"Top features — {m} (mean standardized coefficient over folds; "
            "sign stability = share of folds agreeing in sign):",
            "",
            _fmt(coef.head(8).reset_index(names="feature")),
            "",
        ]

    severity = cfg.severity or next(
        (c for c in SEVERITY_CANDIDATES if c in cohort.subjects.columns), None
    )
    if severity:
        target = cohort.subjects.loc[y == 1, severity]
        rows = []
        for m, X in blocks.items():
            r = regress_severity_cv(X.loc[X.index.intersection(target.index)], target, cov,
                                    seed=cfg.seed)
            rows.append({"block": m, **r})
        sev = pd.DataFrame(rows)
        sev.to_csv(out / "severity_regression.csv", index=False)
        summary["severity"] = {"target": severity, "results": rows}
        md += [
            f"### Symptom severity within ASD: `{severity}` (ridge, 5-fold CV)",
            "",
            _fmt(sev),
            "",
        ]

    (out / "report.md").write_text("\n".join(md) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    return summary
