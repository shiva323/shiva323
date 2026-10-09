"""Resting-state fMRI -> functional-connectivity features.

``connectivity_features`` works on any list of ROI time series (one array of
shape (time, n_rois) per subject), e.g. from fMRIPrep + nilearn maskers or
the preprocessed ABIDE release. Network-level summaries keep the feature
count small enough for cohorts of a few hundred children.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from nilearn.connectome import ConnectivityMeasure


def connectivity_features(
    timeseries: list[np.ndarray],
    subject_ids: list[str],
    roi_networks: list[str] | None = None,
    kind: str = "correlation",
) -> pd.DataFrame:
    """Fisher-z functional connectivity.

    Without ``roi_networks`` every ROI pair is a feature (n_rois*(n_rois-1)/2
    columns — reduce it with PCA before modelling). With a network label per
    ROI (e.g. Yeo-7: 'DMN', 'Salience', ...) features are the mean within- and
    between-network connectivity, which is far easier to interpret."""
    conn = ConnectivityMeasure(kind=kind)
    mats = conn.fit_transform(timeseries)
    if kind == "correlation":
        mats = np.arctanh(np.clip(mats, -0.999999, 0.999999))
    n_rois = mats.shape[1]

    if roi_networks is None:
        iu = np.triu_indices(n_rois, 1)
        cols = [f"fc_{i}_{j}" for i, j in zip(*iu)]
        return pd.DataFrame(mats[:, iu[0], iu[1]], index=subject_ids, columns=cols)

    labels = np.asarray(roi_networks)
    nets = list(dict.fromkeys(roi_networks))
    out = {}
    for a, na in enumerate(nets):
        for nb in nets[a:]:
            ia, ib = np.flatnonzero(labels == na), np.flatnonzero(labels == nb)
            block = mats[:, ia][:, :, ib]
            if na == nb:
                iu = np.triu_indices(len(ia), 1)
                vals = block[:, iu[0], iu[1]].mean(axis=1)
            else:
                vals = block.mean(axis=(1, 2))
            out[f"fc_{na}__{nb}"] = vals
    return pd.DataFrame(out, index=subject_ids)


def fetch_abide_features(data_dir: str | None = None, atlas: str = "rois_cc200"):
    """Download ABIDE I preprocessed ROI time series (C-PAC pipeline) with
    nilearn and return (features, phenotypes). Needs internet (~1-2 GB).

    ABIDE has no microbiome data: use it to build and validate the *brain*
    branch of the model, not for subject-level gut–brain fusion."""
    from nilearn import datasets

    abide = datasets.fetch_abide_pcp(
        data_dir=data_dir,
        derivatives=[atlas],
        pipeline="cpac",
        band_pass_filtering=True,
        global_signal_regression=False,
        quality_checked=True,
    )
    pheno = pd.DataFrame(abide.phenotypic)
    ids = pheno["SUB_ID"].astype(str).tolist()
    feats = connectivity_features(list(abide[atlas]), ids)
    pheno = pheno.set_index(pheno["SUB_ID"].astype(str))
    pheno["group"] = np.where(pheno["DX_GROUP"] == 1, "ASD", "TD")
    return feats, pheno
