"""Report figures (matplotlib, static PNG)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

GROUP_COLORS = {"ASD": "#eb6834", "TD": "#2a78d6"}
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
NEUTRAL = "#c3c2b7"
BAR = "#2a78d6"  # single-series bars
WEIGHT = "#4a3aa7"  # loadings: one hue, sign shown by bar direction (not group colours)

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK_2,
        "axes.titlecolor": INK,
        "axes.titleweight": "bold",
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "font.family": "sans-serif",
    }
)


def _save(fig, path: Path) -> Path:
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def alpha_diversity(alpha: pd.DataFrame, groups: pd.Series, path: Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.0))
    rng = np.random.default_rng(0)
    for ax, metric in zip(axes, ["shannon", "observed"]):
        for i, g in enumerate(GROUP_COLORS):
            vals = alpha.loc[groups == g, metric]
            ax.boxplot(
                vals, positions=[i], widths=0.5, showfliers=False,
                medianprops={"color": INK, "linewidth": 1.5},
                boxprops={"color": AXIS}, whiskerprops={"color": AXIS}, capprops={"color": AXIS},
            )
            ax.scatter(
                i + rng.uniform(-0.15, 0.15, len(vals)), vals, s=10,
                color=GROUP_COLORS[g], alpha=0.7, edgecolor="none",
            )
        ax.set_xticks([0, 1], list(GROUP_COLORS))
        ax.set_title(f"Alpha diversity: {metric}")
        ax.grid(axis="x", visible=False)
    return _save(fig, path)


def ordination(coords: np.ndarray, var: np.ndarray, groups: pd.Series, title: str, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    for g, c in GROUP_COLORS.items():
        m = (groups == g).to_numpy()
        ax.scatter(coords[m, 0], coords[m, 1], s=16, color=c, alpha=0.75,
                   edgecolor=SURFACE, linewidth=0.6, label=g)
    ax.set_xlabel(f"PCo1 ({var[0]:.1%})")
    ax.set_ylabel(f"PCo2 ({var[1]:.1%})")
    ax.set_title(title)
    ax.legend(loc="best")
    return _save(fig, path)


def volcano(res: pd.DataFrame, title: str, path: Path, q: float = 0.05, n_labels: int = 8) -> Path:
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    y = -np.log10(res["q_value"].clip(lower=1e-12))
    sig = res["q_value"] < q
    ax.scatter(res.loc[~sig, "effect_d"], y[~sig], s=10, color=NEUTRAL, edgecolor="none",
               label=f"q ≥ {q}")
    up, down = sig & (res["effect_d"] > 0), sig & (res["effect_d"] < 0)
    ax.scatter(res.loc[up, "effect_d"], y[up], s=18, color=GROUP_COLORS["ASD"],
               edgecolor=SURFACE, linewidth=0.6, label="higher in ASD")
    ax.scatter(res.loc[down, "effect_d"], y[down], s=18, color=GROUP_COLORS["TD"],
               edgecolor=SURFACE, linewidth=0.6, label="higher in TD")
    placed: list[tuple[float, float]] = []
    span_x = float(res["effect_d"].abs().max()) or 1.0
    span_y = float(y.max()) or 1.0
    for _, row in res.head(n_labels).iterrows():
        px, py = row["effect_d"], -np.log10(max(row["q_value"], 1e-12))
        # stack labels of near-coincident points instead of overprinting them
        clash = sum(abs(px - qx) / span_x < 0.15 and abs(py - qy) / span_y < 0.05 for qx, qy in placed)
        placed.append((px, py))
        ax.annotate(row["feature"], (px, py), fontsize=7, color=INK_2,
                    xytext=(5, 2 - 10 * clash), textcoords="offset points",
                    ha="left" if px <= 0 else "right")
    ax.axhline(-np.log10(q), color=AXIS, linewidth=1, linestyle="--")
    ax.axvline(0, color=AXIS, linewidth=1)
    ax.set_xlabel("Adjusted effect (SD units, ASD − TD)")
    ax.set_ylabel("−log10 FDR q")
    ax.set_title(title)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)
    return _save(fig, path)


def scca_mode(result, groups: pd.Series, path: Path, top: int = 8) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), gridspec_kw={"width_ratios": [1.1, 1, 1]})
    ax = axes[0]
    xs, ys = result.x_scores["mode1"], result.y_scores["mode1"]
    for g, c in GROUP_COLORS.items():
        m = (groups.loc[xs.index] == g).to_numpy()
        ax.scatter(xs[m], ys[m], s=14, color=c, alpha=0.75, edgecolor=SURFACE, linewidth=0.6, label=g)
    ax.set_xlabel("Gut score (mode 1)")
    ax.set_ylabel("Brain score (mode 1)")
    ax.set_title(f"sCCA mode 1: r = {result.in_sample_r[0]:.2f}, CV r = {result.cv_r:.2f}")
    ax.legend(loc="best")
    for ax, w, label in [(axes[1], result.x_weights, "Gut weights"), (axes[2], result.y_weights, "Brain weights")]:
        s = w["mode1"][w["mode1"] != 0]
        s = s.reindex(s.abs().sort_values(ascending=False).index)[:top][::-1]
        ax.barh(range(len(s)), s.to_numpy(), color=WEIGHT, height=0.6)
        ax.set_yticks(range(len(s)), s.index, fontsize=7.5, color=INK_2)
        ax.axvline(0, color=AXIS, linewidth=1)
        ax.set_title(f"{label} (non-zero: {(w['mode1'] != 0).sum()})")
        ax.grid(axis="y", visible=False)
    return _save(fig, path)


def model_auc(summary: pd.DataFrame, path: Path) -> Path:
    s = summary.sort_values("AUC_mean")
    fig, ax = plt.subplots(figsize=(5.4, 0.45 * len(s) + 1.2))
    ax.barh(range(len(s)), s["AUC_mean"], xerr=s["AUC_sd"], color=BAR, height=0.55,
            error_kw={"ecolor": INK_2, "elinewidth": 1, "capsize": 2})
    for i, (name, row) in enumerate(s.iterrows()):
        ax.text(row["AUC_mean"] + row["AUC_sd"] + 0.01, i, f"{row['AUC_mean']:.2f}",
                va="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(s)), [f"{n} (n={int(r)})" for n, r in s["n_subjects"].items()],
                  fontsize=8, color=INK_2)
    ax.axvline(0.5, color=AXIS, linestyle="--", linewidth=1)
    ax.set_xlim(0.4, 1.05)
    ax.set_xlabel("ROC-AUC, nested CV (mean ± SD over repeats)")
    ax.set_title("ASD vs TD classification by modality")
    ax.grid(axis="y", visible=False)
    return _save(fig, path)


def mediation_diagram(med: dict, x: str, m: str, y: str, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6.4, 2.8))
    ax.set_axis_off()
    ax.grid(False)
    nodes = {x: (0.14, 0.25), m: (0.5, 0.8), y: (0.86, 0.25)}
    boxes = {}
    for name, (px, py) in nodes.items():
        boxes[name] = ax.text(px, py, name, ha="center", va="center", fontsize=9, color=INK,
                              bbox={"boxstyle": "round,pad=0.5", "facecolor": SURFACE,
                                    "edgecolor": AXIS})

    def arrow(a, b, label, offset):
        (x0, y0), (x1, y1) = nodes[a], nodes[b]
        # patchA/patchB clip the arrow at the label boxes, whatever their width
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops={"arrowstyle": "-|>", "color": INK_2, "lw": 1.2,
                                "patchA": boxes[a].get_bbox_patch(),
                                "patchB": boxes[b].get_bbox_patch(),
                                "shrinkA": 4, "shrinkB": 4})
        ax.text((x0 + x1) / 2 + offset[0], (y0 + y1) / 2 + offset[1], label,
                ha="center", fontsize=8.5, color=INK_2)

    arrow(x, m, f"a = {med['a_path']:.2f}", (-0.08, 0.02))
    arrow(m, y, f"b = {med['b_path']:.2f}", (0.08, 0.02))
    arrow(x, y, f"c' = {med['direct_c_prime']:.2f}", (0, -0.12))
    ax.set_title(
        f"Indirect a×b = {med['indirect_ab']:.3f} "
        f"[95% CI {med['indirect_ci_low']:.3f}, {med['indirect_ci_high']:.3f}]",
        fontsize=10,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return _save(fig, path)
