"""Covariate-adjusted univariate tests and mediation analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def _ols(design: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """OLS coefficients, standard errors and residual degrees of freedom."""
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    dof = design.shape[0] - np.linalg.matrix_rank(design)
    sigma2 = resid @ resid / dof
    cov = sigma2 * np.linalg.pinv(design.T @ design)
    return beta, np.sqrt(np.diag(cov)), dof


def group_differences(
    features: pd.DataFrame,
    y: pd.Series,
    covariates: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """For every feature fit ``feature ~ group + covariates`` and test the
    group coefficient. Effect size is the adjusted difference in SD units
    (Cohen's d analogue); q-values are Benjamini–Hochberg FDR."""
    idx = features.index
    g = y.loc[idx].to_numpy(float)
    cols = [np.ones(len(idx)), g]
    if covariates is not None and covariates.shape[1]:
        cols.append(covariates.loc[idx].to_numpy(float))
    design = np.column_stack(cols)

    X = features.to_numpy(float)
    sd = X.std(axis=0, ddof=1)
    rows = []
    for j, name in enumerate(features.columns):
        if sd[j] == 0:
            continue
        beta, se, dof = _ols(design, X[:, j])
        t = beta[1] / se[1]
        rows.append(
            {
                "feature": name,
                "effect_d": beta[1] / sd[j],
                "t": t,
                "p_value": 2 * stats.t.sf(abs(t), dof),
                "mean_ASD": X[g == 1, j].mean(),
                "mean_TD": X[g == 0, j].mean(),
            }
        )
    out = pd.DataFrame(rows)
    out["q_value"] = stats.false_discovery_control(out["p_value"].to_numpy())
    return out.sort_values("p_value").reset_index(drop=True)


def correlate_blocks(
    A: pd.DataFrame, B: pd.DataFrame, method: str = "spearman"
) -> pd.DataFrame:
    """All pairwise correlations between two feature blocks, FDR-corrected.
    Use residualized blocks so covariates do not drive the associations."""
    idx = A.index.intersection(B.index)
    a, b = A.loc[idx], B.loc[idx]
    fn = stats.spearmanr if method == "spearman" else stats.pearsonr
    rows = []
    for ca in a.columns:
        for cb in b.columns:
            r, p = fn(a[ca], b[cb])
            rows.append({"feature_a": ca, "feature_b": cb, "r": r, "p_value": p})
    out = pd.DataFrame(rows)
    out["q_value"] = stats.false_discovery_control(out["p_value"].to_numpy())
    return out.sort_values("p_value").reset_index(drop=True)


def mediation(
    x: pd.Series,
    m: pd.Series,
    y: pd.Series,
    covariates: pd.DataFrame | None = None,
    n_boot: int = 2000,
    seed: int = 0,
) -> dict[str, float]:
    """Single-mediator model X -> M -> Y with bootstrap CI of the indirect effect.

    a-path:  M ~ X + C
    b-path:  Y ~ M + X + C   (c' = direct effect of X)
    indirect = a * b, total = c' + a * b

    Variables are z-scored so effects are in SD units. A cross-sectional
    mediation is consistent with — not proof of — a causal pathway.
    """
    idx = x.index.intersection(m.index).intersection(y.index)
    z = lambda s: ((s - s.mean()) / s.std(ddof=1)).to_numpy(float)  # noqa: E731
    xv, mv, yv = z(x.loc[idx]), z(m.loc[idx]), z(y.loc[idx])
    C = (
        covariates.loc[idx].to_numpy(float)
        if covariates is not None and covariates.shape[1]
        else np.empty((len(idx), 0))
    )
    ones = np.ones(len(idx))

    def paths(rows: np.ndarray) -> tuple[float, float, float]:
        da = np.column_stack([ones[rows], xv[rows], C[rows]])
        a = np.linalg.lstsq(da, mv[rows], rcond=None)[0][1]
        db = np.column_stack([ones[rows], mv[rows], xv[rows], C[rows]])
        coef = np.linalg.lstsq(db, yv[rows], rcond=None)[0]
        return a, coef[1], coef[2]

    all_rows = np.arange(len(idx))
    a, b, c_prime = paths(all_rows)
    rng = np.random.default_rng(seed)
    boot = np.array(
        [np.prod(paths(rng.choice(all_rows, len(all_rows)))[:2]) for _ in range(n_boot)]
    )
    lo, hi = np.percentile(boot, [2.5, 97.5])
    # Two-sided bootstrap p-value for the indirect effect.
    p = 2 * min((boot <= 0).mean(), (boot >= 0).mean())
    total = c_prime + a * b
    return {
        "n": len(idx),
        "a_path": a,
        "b_path": b,
        "direct_c_prime": c_prime,
        "indirect_ab": a * b,
        "indirect_ci_low": lo,
        "indirect_ci_high": hi,
        "indirect_p_boot": max(p, 1 / n_boot),
        "proportion_mediated": (a * b) / total if abs(total) > 1e-8 else np.nan,
    }
