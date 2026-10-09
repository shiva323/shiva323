"""Predictive models evaluated with leakage-free (nested) cross-validation.

Three strategies are compared on the same outer folds:

* single modality   one classifier per block (microbiome, brain, metabolites)
* early fusion      all blocks concatenated (only subjects with every block)
* late fusion       one classifier per block, then a logistic meta-learner on
                    their out-of-fold logits (stacking). A subject missing a
                    block contributes a neutral logit of 0 for it, so children
                    without EEG / MRI still get a prediction.

Covariates (age, sex, diet, site, ...) are regressed out *inside* each
training fold, never on the full data set.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler


class ConfoundRegressor(BaseEstimator, TransformerMixin):
    """Expects the last ``n_confounds`` columns of X to be confounds; fits
    feature ~ confounds on the training data and returns residual features."""

    def __init__(self, n_confounds: int = 0):
        self.n_confounds = n_confounds

    def _split(self, X):
        X = np.asarray(X, float)
        if self.n_confounds == 0:
            return X, np.ones((len(X), 1))
        feats, conf = X[:, : -self.n_confounds], X[:, -self.n_confounds :]
        return feats, np.column_stack([np.ones(len(X)), conf])

    def fit(self, X, y=None):
        feats, conf = self._split(X)
        self.beta_ = np.linalg.lstsq(conf, feats, rcond=None)[0]
        return self

    def transform(self, X):
        feats, conf = self._split(X)
        return feats - conf @ self.beta_


def _with_confounds(X: pd.DataFrame, conf: pd.DataFrame | None) -> np.ndarray:
    if conf is None or conf.shape[1] == 0:
        return X.to_numpy(float)
    return np.column_stack([X.to_numpy(float), conf.loc[X.index].to_numpy(float)])


def make_classifier(n_confounds: int, block_sizes: list[int] | None = None) -> GridSearchCV:
    """Covariate regression -> z-score -> L2 logistic regression, with the
    regularisation strength tuned in an inner 5-fold CV (the nested loop).

    ``block_sizes`` (early fusion) down-weights each block by 1/sqrt(size) so
    the block with most features does not dominate the concatenation."""
    steps = [("confounds", ConfoundRegressor(n_confounds)), ("scale", StandardScaler())]
    if block_sizes:
        w = np.concatenate([np.full(s, 1 / np.sqrt(s)) for s in block_sizes])
        steps.append(("block_weight", FunctionTransformer(lambda Z, w=w: Z * w)))
    steps.append(("clf", LogisticRegression(max_iter=5000, class_weight="balanced")))
    return GridSearchCV(
        Pipeline(steps),
        {"clf__C": np.logspace(-3, 1, 9)},
        cv=StratifiedKFold(5, shuffle=True, random_state=0),
        scoring="roc_auc",
    )


def _logit(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return np.log(p / (1 - p))


def _zscore(x: np.ndarray, ref: np.ndarray | None = None) -> np.ndarray:
    ref = x if ref is None else ref
    return (x - ref.mean()) / (ref.std() or 1.0)


@dataclass
class FusionResult:
    auc: pd.DataFrame            # rows = repeats, columns = model
    auc_complete: pd.DataFrame   # same, restricted to subjects with every block
    balanced_accuracy: pd.DataFrame
    oof_proba: pd.DataFrame      # subject x model, averaged over repeats
    coefficients: dict[str, pd.DataFrame]  # per block: mean coef, sign stability
    meta_weights: pd.DataFrame   # late-fusion weight per block (z-scored logits) per fold
    n_subjects: dict[str, int]

    def summary(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "AUC_mean": self.auc.mean(),
                "AUC_sd": self.auc.std(ddof=1) if len(self.auc) > 1 else 0.0,
                "AUC_complete_cases": self.auc_complete.mean(),
                "balanced_acc": self.balanced_accuracy.mean(),
                "n_subjects": pd.Series(self.n_subjects),
            }
        ).sort_values("AUC_mean", ascending=False)


def run_fusion_cv(
    blocks: dict[str, pd.DataFrame],
    y: pd.Series,
    confounds: pd.DataFrame | None = None,
    n_splits: int = 5,
    n_repeats: int = 3,
    seed: int = 0,
) -> FusionResult:
    subjects = y.index
    n_conf = 0 if confounds is None else confounds.shape[1]
    complete = subjects
    for X in blocks.values():
        complete = complete.intersection(X.index)
    early_X = pd.concat([X.loc[complete] for X in blocks.values()], axis=1)
    early_sizes = [X.shape[1] for X in blocks.values()]

    names = list(blocks) + (["early_fusion"] if len(blocks) > 1 else []) + ["late_fusion"]
    auc_rows, auc_cc_rows, bacc_rows, proba_runs, meta_rows = [], [], [], [], []
    coefs: dict[str, list[np.ndarray]] = {m: [] for m in blocks}

    for rep in range(n_repeats):
        outer = StratifiedKFold(n_splits, shuffle=True, random_state=seed + rep)
        proba = pd.DataFrame(np.nan, index=subjects, columns=names)
        for tr, te in outer.split(np.zeros(len(subjects)), y):
            s_tr, s_te = subjects[tr], subjects[te]
            meta_tr = np.zeros((len(s_tr), len(blocks)))
            meta_te = np.zeros((len(s_te), len(blocks)))
            for k, (m, X) in enumerate(blocks.items()):
                a_tr, a_te = s_tr.intersection(X.index), s_te.intersection(X.index)
                Xtr = _with_confounds(X.loc[a_tr], confounds)
                model = make_classifier(n_conf)
                inner = StratifiedKFold(5, shuffle=True, random_state=seed + rep)
                p_inner = cross_val_predict(
                    clone(model), Xtr, y.loc[a_tr], cv=inner, method="predict_proba"
                )[:, 1]
                model.fit(Xtr, y.loc[a_tr])
                coefs[m].append(model.best_estimator_.named_steps["clf"].coef_[0])
                # Base models differ wildly in calibration (logit SD 0.3 vs 3) and
                # the refit model may pick another C than the inner fits, so each
                # set of logits is z-scored against the model that produced it.
                # A missing block then sits at 0 = average, uninformative.
                meta_tr[s_tr.get_indexer(a_tr), k] = _zscore(_logit(p_inner))
                if len(a_te):
                    ref = _logit(model.predict_proba(Xtr)[:, 1])
                    p_te = model.predict_proba(_with_confounds(X.loc[a_te], confounds))[:, 1]
                    proba.loc[a_te, m] = p_te
                    meta_te[s_te.get_indexer(a_te), k] = _zscore(_logit(p_te), ref)
            meta = LogisticRegression(C=1.0, class_weight="balanced").fit(meta_tr, y.loc[s_tr])
            meta_rows.append(dict(zip(blocks, meta.coef_[0])))
            proba.loc[s_te, "late_fusion"] = meta.predict_proba(meta_te)[:, 1]

            if len(blocks) > 1:
                e_tr, e_te = s_tr.intersection(complete), s_te.intersection(complete)
                model = make_classifier(n_conf, early_sizes)
                model.fit(_with_confounds(early_X.loc[e_tr], confounds), y.loc[e_tr])
                if len(e_te):
                    proba.loc[e_te, "early_fusion"] = model.predict_proba(
                        _with_confounds(early_X.loc[e_te], confounds)
                    )[:, 1]

        auc, auc_cc, bacc = {}, {}, {}
        for m in names:
            ok = proba[m].notna()
            auc[m] = roc_auc_score(y[ok], proba.loc[ok, m])
            auc_cc[m] = roc_auc_score(y.loc[complete], proba.loc[complete, m])
            bacc[m] = balanced_accuracy_score(y[ok], proba.loc[ok, m] >= 0.5)
        auc_rows.append(auc)
        auc_cc_rows.append(auc_cc)
        bacc_rows.append(bacc)
        proba_runs.append(proba)

    coef_tables = {}
    for m, X in blocks.items():
        C = np.array(coefs[m])
        mean = C.mean(axis=0)
        coef_tables[m] = (
            pd.DataFrame(
                {
                    "mean_coef": mean,
                    "sign_stability": (np.sign(C) == np.sign(mean)).mean(axis=0),
                },
                index=X.columns,
            )
            .assign(abs_coef=lambda d: d["mean_coef"].abs())
            .sort_values("abs_coef", ascending=False)
            .drop(columns="abs_coef")
        )

    n_subj = {m: len(X.index.intersection(subjects)) for m, X in blocks.items()}
    n_subj["late_fusion"] = len(subjects)
    if len(blocks) > 1:
        n_subj["early_fusion"] = len(complete)
    return FusionResult(
        auc=pd.DataFrame(auc_rows),
        auc_complete=pd.DataFrame(auc_cc_rows),
        balanced_accuracy=pd.DataFrame(bacc_rows),
        oof_proba=sum(proba_runs) / len(proba_runs),
        coefficients=coef_tables,
        meta_weights=pd.DataFrame(meta_rows),
        n_subjects=n_subj,
    )


def regress_severity_cv(
    X: pd.DataFrame,
    target: pd.Series,
    confounds: pd.DataFrame | None = None,
    n_splits: int = 5,
    seed: int = 0,
) -> dict[str, float]:
    """Predict a continuous score (e.g. SRS / ADOS CSS / CARS) with ridge
    regression; returns out-of-fold Pearson r, R^2 and MAE."""
    idx = X.index.intersection(target.dropna().index)
    n_conf = 0 if confounds is None else confounds.shape[1]
    model = Pipeline(
        [
            ("confounds", ConfoundRegressor(n_conf)),
            ("scale", StandardScaler()),
            ("ridge", RidgeCV(alphas=np.logspace(-2, 4, 25))),
        ]
    )
    yv = target.loc[idx].to_numpy(float)
    pred = cross_val_predict(
        model,
        _with_confounds(X.loc[idx], confounds),
        yv,
        cv=KFold(n_splits, shuffle=True, random_state=seed),
    )
    ss_res = ((yv - pred) ** 2).sum()
    ss_tot = ((yv - yv.mean()) ** 2).sum()
    return {
        "n": len(idx),
        "r": float(np.corrcoef(yv, pred)[0, 1]),
        "R2": float(1 - ss_res / ss_tot),
        "MAE": float(np.abs(yv - pred).mean()),
    }
