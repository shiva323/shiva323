"""Smoke and recovery tests: each method must find what the simulator planted."""

import sys
from pathlib import Path

import mne
import numpy as np
import pandas as pd
import pytest
from scipy.signal import detrend

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gutbrain import eeg, fmri, microbiome as mb  # noqa: E402
from gutbrain.data import load_cohort, save_cohort  # noqa: E402
from gutbrain.integration import sparse_cca  # noqa: E402
from gutbrain.pipeline import PipelineConfig, run_pipeline  # noqa: E402
from gutbrain.simulate import simulate_cohort  # noqa: E402
from gutbrain.stats import group_differences, mediation  # noqa: E402


@pytest.fixture(scope="module")
def cohort():
    return simulate_cohort(n_asd=60, n_td=60, seed=3)


def test_clr_and_permanova(cohort):
    clr = mb.clr(mb.filter_taxa(cohort.microbiome))
    assert np.allclose(clr.sum(axis=1), 0)
    perm = mb.permanova(mb.aitchison_distance(clr), cohort.subjects["group"].to_numpy(), n_perm=199)
    assert perm["p_value"] < 0.05


def test_differential_abundance_finds_planted_taxa(cohort):
    clr = mb.clr(mb.filter_taxa(cohort.microbiome))
    res = group_differences(clr, cohort.y).set_index("feature")
    assert res.loc["Bifidobacterium", "effect_d"] < 0
    assert res.loc["Desulfovibrio", "effect_d"] > 0
    assert res.loc["Bifidobacterium", "q_value"] < 0.05


def test_sparse_cca_recovers_gut_brain_mode(cohort):
    both = cohort.brain.index
    clr = mb.clr(mb.filter_taxa(cohort.microbiome)).loc[both]
    res = sparse_cca(clr, cohort.brain, n_components=1, n_perm=50)
    top_gut = res.x_weights["mode1"].abs().sort_values(ascending=False).index[:3]
    top_brain = res.y_weights["mode1"].abs().sort_values(ascending=False).index[:2]
    assert {"Bifidobacterium", "Desulfovibrio"} & set(top_gut)
    assert {"wpli_alpha_frontal", "aperiodic_exponent"} & set(top_brain)
    assert res.perm_p < 0.05 and res.cv_r > 0.3


def test_mediation_recovers_indirect_path(cohort):
    clr = mb.clr(cohort.microbiome)
    out = mediation(
        clr["Bifidobacterium"],
        cohort.brain["wpli_alpha_frontal"],
        cohort.subjects["SRS_total"],
        pd.DataFrame({"ASD": cohort.y}),
        n_boot=500,
    )
    assert out["indirect_ci_high"] < 0  # more Bifidobacterium -> higher wPLI -> lower SRS


def _synthetic_raw(alpha_hz=10.0, seconds=60, sfreq=250.0, seed=0):
    rng = np.random.default_rng(seed)
    name = "colin27_1020" if "colin27_1020" in mne.channels.get_builtin_montages() else "standard_1020"
    montage = mne.channels.make_standard_montage(name)
    chs = ["Fp1", "Fp2", "F3", "F4", "Fz", "C3", "C4", "Cz", "T7", "T8", "P3", "P4", "Pz", "O1", "O2"]
    t = np.arange(int(seconds * sfreq)) / sfreq
    brown = detrend(np.cumsum(rng.normal(size=(len(chs), len(t))), axis=1), axis=1)
    data = brown / brown.std() + rng.normal(size=brown.shape) * 0.5
    posterior = [chs.index(c) for c in ("P3", "P4", "Pz", "O1", "O2")]
    data[posterior] += 3 * np.sin(2 * np.pi * alpha_hz * t)
    info = mne.create_info(chs, sfreq, "eeg")
    raw = mne.io.RawArray(data * 10e-6, info, verbose="error")
    raw.set_montage(montage)
    return raw


def test_eeg_features():
    feats = eeg.extract_features(_synthetic_raw(alpha_hz=10.0), line_freq=50.0)
    assert abs(feats["alpha_peak_freq"] - 10.0) <= 0.5
    assert feats["relpow_alpha_occipital"] > feats["relpow_alpha_frontal"]
    assert 0 <= feats["wpli_alpha_global"] <= 1
    assert feats["aperiodic_exponent"] > 0


def test_region_mapping_by_position():
    info = mne.create_info([f"E{i}" for i in range(1, 33)], 250.0, "eeg")
    montage = mne.channels.make_standard_montage("GSN-HydroCel-32")
    info.set_montage(mne.channels.make_dig_montage(
        ch_pos=dict(zip(info.ch_names, montage.get_positions()["ch_pos"].values())),
        coord_frame="head",
    ))
    regions = eeg.channel_regions(info)
    assert set(regions.values()) <= set(eeg.REGIONS)
    assert len(set(regions.values())) >= 4


def test_fmri_network_features():
    rng = np.random.default_rng(0)
    ts = [rng.normal(size=(120, 6)) for _ in range(4)]
    feats = fmri.connectivity_features(ts, ["a", "b", "c", "d"], ["DMN"] * 3 + ["SAL"] * 3)
    assert list(feats.columns) == ["fc_DMN__DMN", "fc_DMN__SAL", "fc_SAL__SAL"]


def test_end_to_end(tmp_path, cohort):
    save_cohort(cohort, tmp_path / "data")
    loaded = load_cohort(tmp_path / "data")
    cfg = PipelineConfig(
        covariates=["age", "sex", "site", "fiber_score"],
        mediation=[("Bifidobacterium", "wpli_alpha_frontal", "SRS_total")],
        n_perm=20,
        cv_repeats=1,
    )
    summary = run_pipeline(loaded, tmp_path / "out", cfg)
    assert (tmp_path / "out" / "report.md").exists()
    assert summary["classification"]["late_fusion"]["AUC_mean"] > 0.5
    assert summary["scca"]["perm_p"] < 0.1
