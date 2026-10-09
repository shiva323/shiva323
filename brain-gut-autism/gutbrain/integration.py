"""Sparse canonical correlation analysis (sCCA) between two feature blocks.

Implements the penalized matrix decomposition of Witten, Tibshirani & Hastie
(2009, Biostatistics). With hundreds of taxa and dozens of brain features the
classic CCA overfits perfectly; L1 penalties select a small set of taxa and a
small set of brain features whose weighted sums co-vary across children — a
data-driven "gut–brain axis mode".
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold


def _soft(x: np.ndarray, delta: float) -> np.ndarray:
    return np.sign(x) * np.maximum(np.abs(x) - delta, 0.0)


def _l1_project(a: np.ndarray, c: float) -> np.ndarray:
    """Unit-L2 vector S(a, delta)/||S(a, delta)|| with L1 norm <= c."""
    norm = np.linalg.norm(a)
    if norm == 0:
        return a
    u = a / norm
    if np.abs(u).sum() <= c:
        return u
    lo, hi = 0.0, np.abs(a).max()
    for _ in range(60):  # binary search for the threshold
        mid = (lo + hi) / 2
        s = _soft(a, mid)
        s_norm = np.linalg.norm(s)
        if s_norm == 0 or np.abs(s).sum() / s_norm < c:
            hi = mid
        else:
            lo = mid
    s = _soft(a, lo)
    return s / np.linalg.norm(s)


def _standardize(train: np.ndarray, *others: np.ndarray) -> list[np.ndarray]:
    mu, sd = train.mean(axis=0), train.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    return [(m - mu) / sd for m in (train, *others)]


def _pmd(X, Y, cx, cy, n_components=1, n_iter=200, tol=1e-7):
    K = X.T @ Y
    U, V = [], []
    for _ in range(n_components):
        v = np.linalg.svd(K, full_matrices=False)[2][0]
        u = np.zeros(K.shape[0])
        for _ in range(n_iter):
            u_new = _l1_project(K @ v, cx)
            v_new = _l1_project(K.T @ u_new, cy)
            done = np.abs(u_new - u).sum() + np.abs(v_new - v).sum() < tol
            u, v = u_new, v_new
            if done:
                break
        d = u @ K @ v
        K = K - d * np.outer(u, v)  # deflate before the next component
        U.append(u)
        V.append(v)
    return np.array(U).T, np.array(V).T


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


@dataclass
class SCCAResult:
    x_weights: pd.DataFrame
    y_weights: pd.DataFrame
    x_scores: pd.DataFrame
    y_scores: pd.DataFrame
    in_sample_r: list[float]
    cv_r: float
    perm_p: float
    penalty_x: float
    penalty_y: float


def sparse_cca(
    X: pd.DataFrame,
    Y: pd.DataFrame,
    n_components: int = 2,
    penalty_grid: tuple[float, ...] = (0.1, 0.2, 0.3, 0.5, 0.7),
    n_folds: int = 5,
    n_perm: int = 500,
    seed: int = 0,
) -> SCCAResult:
    """Fit sCCA with penalties chosen by cross-validated canonical correlation.

    ``penalty_grid`` values are fractions of sqrt(p): small = sparser.
    ``cv_r`` (out-of-sample correlation of the first mode) is the honest effect
    size; ``perm_p`` tests the first in-sample correlation against row
    permutations of Y with the selected penalties.
    Residualize both blocks on covariates (age, sex, diet, site) beforehand.
    """
    idx = X.index.intersection(Y.index)
    Xv, Yv = X.loc[idx].to_numpy(float), Y.loc[idx].to_numpy(float)
    p, q = Xv.shape[1], Yv.shape[1]
    to_c = lambda f, dim: max(1.0, f * np.sqrt(dim))  # noqa: E731
    folds = list(KFold(n_folds, shuffle=True, random_state=seed).split(Xv))

    def cv_score(fx: float, fy: float) -> float:
        rs = []
        for tr, te in folds:
            Xtr, Xte = _standardize(Xv[tr], Xv[te])
            Ytr, Yte = _standardize(Yv[tr], Yv[te])
            u, v = _pmd(Xtr, Ytr, to_c(fx, p), to_c(fy, q))
            rs.append(_corr(Xte @ u[:, 0], Yte @ v[:, 0]))
        return float(np.mean(rs))

    grid = [(fx, fy) for fx in penalty_grid for fy in penalty_grid]
    scores = [cv_score(fx, fy) for fx, fy in grid]
    fx, fy = grid[int(np.argmax(scores))]
    cx, cy = to_c(fx, p), to_c(fy, q)

    Xs, Ys = _standardize(Xv)[0], _standardize(Yv)[0]
    U, V = _pmd(Xs, Ys, cx, cy, n_components=n_components)
    xs, ys = Xs @ U, Ys @ V
    in_r = [_corr(xs[:, k], ys[:, k]) for k in range(n_components)]

    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_perm):
        Yp = Ys[rng.permutation(len(Ys))]
        u, v = _pmd(Xs, Yp, cx, cy)
        null.append(_corr(Xs @ u[:, 0], Yp @ v[:, 0]))
    perm_p = (1 + np.sum(np.array(null) >= in_r[0])) / (n_perm + 1)

    comp = [f"mode{k + 1}" for k in range(n_components)]
    return SCCAResult(
        x_weights=pd.DataFrame(U, index=X.columns, columns=comp),
        y_weights=pd.DataFrame(V, index=Y.columns, columns=comp),
        x_scores=pd.DataFrame(xs, index=idx, columns=comp),
        y_scores=pd.DataFrame(ys, index=idx, columns=comp),
        in_sample_r=in_r,
        cv_r=max(scores),
        perm_p=float(perm_p),
        penalty_x=fx,
        penalty_y=fy,
    )
