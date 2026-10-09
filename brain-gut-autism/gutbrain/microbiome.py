"""Microbiome preprocessing and ecology statistics.

Sequencing counts are *compositional*: only ratios between taxa carry
information, so features are analysed after a centred log-ratio (CLR)
transform and beta diversity uses the Aitchison distance (Euclidean on CLR).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform


def filter_taxa(
    counts: pd.DataFrame, min_prevalence: float = 0.1, min_rel_abundance: float = 1e-4
) -> pd.DataFrame:
    """Keep taxa present in >= ``min_prevalence`` of samples with a mean
    relative abundance >= ``min_rel_abundance``. Rare taxa are noisy and
    inflate the multiple-testing burden."""
    rel = counts.div(counts.sum(axis=1), axis=0)
    keep = ((counts > 0).mean() >= min_prevalence) & (rel.mean() >= min_rel_abundance)
    return counts.loc[:, keep]


def clr(counts: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    """Centred log-ratio transform; the pseudocount handles zero counts."""
    logx = np.log(counts + pseudocount)
    return logx.sub(logx.mean(axis=1), axis=0)


def alpha_diversity(counts: pd.DataFrame) -> pd.DataFrame:
    """Observed richness, Shannon and inverse Simpson per sample.

    Richness depends on sequencing depth; rarefy or include depth as a
    covariate when depths differ a lot between groups."""
    rel = counts.div(counts.sum(axis=1), axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        shannon = -(rel * np.log(rel)).fillna(0).sum(axis=1)
    return pd.DataFrame(
        {
            "observed": (counts > 0).sum(axis=1),
            "shannon": shannon,
            "inv_simpson": 1.0 / (rel**2).sum(axis=1),
            "depth": counts.sum(axis=1),
        }
    )


def residualize(X: pd.DataFrame, covariates: pd.DataFrame | None) -> pd.DataFrame:
    """Regress covariates out of every column of ``X`` (OLS with intercept)."""
    if covariates is None or covariates.shape[1] == 0:
        return X - X.mean()
    C = np.column_stack([np.ones(len(X)), covariates.loc[X.index].to_numpy(float)])
    beta, *_ = np.linalg.lstsq(C, X.to_numpy(float), rcond=None)
    return pd.DataFrame(X.to_numpy(float) - C @ beta, index=X.index, columns=X.columns)


def permanova(
    distance: np.ndarray, groups: np.ndarray, n_perm: int = 999, seed: int = 0
) -> dict[str, float]:
    """One-way PERMANOVA (Anderson 2001): pseudo-F and permutation p-value."""
    rng = np.random.default_rng(seed)
    d2 = np.asarray(distance) ** 2
    n = len(groups)
    labels = np.unique(groups)
    ss_total = d2[np.triu_indices(n, 1)].sum() / n

    def ss_within(g: np.ndarray) -> float:
        total = 0.0
        for lab in labels:
            idx = np.flatnonzero(g == lab)
            sub = d2[np.ix_(idx, idx)]
            total += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
        return total

    def pseudo_f(ssw: float) -> float:
        return ((ss_total - ssw) / (len(labels) - 1)) / (ssw / (n - len(labels)))

    ssw = ss_within(groups)
    observed = pseudo_f(ssw)
    null = np.array([pseudo_f(ss_within(rng.permutation(groups))) for _ in range(n_perm)])
    return {
        "pseudo_F": float(observed),
        "R2": float((ss_total - ssw) / ss_total),
        "p_value": float((1 + (null >= observed).sum()) / (n_perm + 1)),
    }


def aitchison_distance(clr_values: pd.DataFrame) -> np.ndarray:
    return squareform(pdist(clr_values.to_numpy(float), metric="euclidean"))


def pcoa(distance: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Principal coordinates; returns coordinates and explained-variance ratios."""
    n = distance.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (distance**2) @ J
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    pos = vals > 0
    coords = vecs[:, :n_components] * np.sqrt(np.maximum(vals[:n_components], 0))
    return coords, vals[:n_components] / vals[pos].sum()
