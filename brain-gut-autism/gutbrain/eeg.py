"""Resting-state EEG -> subject-level features for the gut–brain models.

Features per subject
--------------------
* relative band power (delta, theta, alpha, beta, gamma) per scalp region
* frontal theta/beta ratio and individual alpha peak frequency (posterior)
* aperiodic exponent and offset (1/f slope, a proxy of E/I balance)
* weighted phase-lag index (wPLI) connectivity: global mean per band and the
  node strength of every region per band (wPLI is insensitive to volume
  conduction at zero lag, which matters for children's sensor-space EEG)

Any format MNE can read works (EDF, BDF, EEGLAB .set, BrainVision, FIF, EGI).
For source-space features, apply the same functions to label time courses
extracted with ``mne.extract_label_time_course`` (e.g. Desikan atlas).
"""

from __future__ import annotations

from pathlib import Path

import mne
import numpy as np
import pandas as pd

BANDS = {
    "delta": (1.0, 4.0),
    "theta": (4.0, 8.0),
    "alpha": (8.0, 13.0),
    "beta": (13.0, 30.0),
    "gamma": (30.0, 45.0),
}
REGIONS = ("frontal", "central", "temporal", "parietal", "occipital")

_PREFIXES = [  # two-letter prefixes are checked before one-letter ones
    ("FP", "frontal"), ("AF", "frontal"), ("FT", "temporal"), ("TP", "temporal"),
    ("FC", "central"), ("CP", "central"), ("PO", "occipital"),
    ("F", "frontal"), ("C", "central"), ("T", "temporal"), ("P", "parietal"),
    ("O", "occipital"), ("I", "occipital"),
]


def load_raw(path: str | Path) -> mne.io.BaseRaw:
    return mne.io.read_raw(path, preload=True, verbose="error")


def preprocess(
    raw: mne.io.BaseRaw,
    l_freq: float = 1.0,
    h_freq: float = 45.0,
    line_freq: float | None = 50.0,
    sfreq: float = 250.0,
) -> mne.io.BaseRaw:
    """Minimal, reproducible cleaning: EEG channels, notch, band-pass,
    resample, average reference. Use ICA (or ASR) beforehand for heavily
    contaminated paediatric recordings."""
    raw = raw.copy().pick("eeg", exclude="bads")
    if line_freq and line_freq < raw.info["sfreq"] / 2:
        raw.notch_filter(line_freq, verbose="error")
    raw.filter(l_freq, h_freq, verbose="error")
    if raw.info["sfreq"] > sfreq:
        raw.resample(sfreq, verbose="error")
    raw.set_eeg_reference("average", verbose="error")
    return raw


def make_epochs(
    raw: mne.io.BaseRaw, length: float = 2.0, reject_uv: float = 150.0
) -> mne.Epochs:
    epochs = mne.make_fixed_length_epochs(raw, duration=length, preload=True, verbose="error")
    epochs.drop_bad(reject={"eeg": reject_uv * 1e-6}, verbose="error")
    if len(epochs) < 10:
        raise ValueError(f"only {len(epochs)} clean epochs left; check the recording")
    return epochs


def channel_regions(info: mne.Info) -> dict[str, str]:
    """Map channels to scalp regions from 10-20 style names, falling back to
    sensor positions (e.g. EGI 'E1'..'E128' nets)."""
    out: dict[str, str] = {}
    for ch in info.ch_names:
        name = ch.upper().replace("EEG", "").strip(" -_")
        region = next((r for p, r in _PREFIXES if name.startswith(p)), None)
        if region is None:
            region = _region_from_position(info, ch)
        if region:
            out[ch] = region
    return out


def _region_from_position(info: mne.Info, ch: str) -> str | None:
    loc = info["chs"][info.ch_names.index(ch)]["loc"][:3]
    if not np.all(np.isfinite(loc)) or np.allclose(loc, 0):
        return None
    pos = np.array([c["loc"][:3] for c in info["chs"]])
    pos = pos[np.isfinite(pos).all(axis=1) & ~np.all(pos == 0, axis=1)]
    center = pos.mean(axis=0)
    radius = np.linalg.norm(pos - center, axis=1).max()
    x, y = (loc[:2] - center[:2]) / radius
    if y > 0.35:
        return "frontal"
    if y < -0.6:
        return "occipital"
    if abs(x) > 0.6:
        return "temporal"
    if y < -0.2:
        return "parietal"
    return "central"


def spectral_features(epochs: mne.Epochs, regions: dict[str, str]) -> dict[str, float]:
    sfreq = epochs.info["sfreq"]
    spec = epochs.compute_psd(
        method="welch", fmin=1.0, fmax=45.0, n_fft=int(2 * sfreq), verbose="error"
    )
    psd, freqs = spec.get_data(return_freqs=True)
    psd = psd.mean(axis=0)  # (channels, freqs), averaged over epochs
    ch_names = epochs.ch_names
    total = psd.sum(axis=1, keepdims=True)

    feats: dict[str, float] = {}
    for band, (lo, hi) in BANDS.items():
        mask = (freqs >= lo) & (freqs < hi)
        rel = psd[:, mask].sum(axis=1) / total[:, 0]
        for region in REGIONS:
            idx = [i for i, ch in enumerate(ch_names) if regions.get(ch) == region]
            if idx:
                feats[f"relpow_{band}_{region}"] = float(rel[idx].mean())

    def band_power(region: str, band: str) -> float:
        idx = [i for i, ch in enumerate(ch_names) if regions.get(ch) == region]
        lo, hi = BANDS[band]
        mask = (freqs >= lo) & (freqs < hi)
        return float(psd[np.ix_(idx, mask)].sum(axis=1).mean()) if idx else np.nan

    feats["theta_beta_ratio_frontal"] = band_power("frontal", "theta") / band_power(
        "frontal", "beta"
    )
    post = [i for i, ch in enumerate(ch_names) if regions.get(ch) in ("parietal", "occipital")]
    post = post or list(range(len(ch_names)))
    amask = (freqs >= 7.0) & (freqs <= 13.0)
    feats["alpha_peak_freq"] = float(freqs[amask][psd[np.ix_(post, amask)].mean(0).argmax()])

    # Aperiodic 1/f fit in log-log space, excluding the alpha peak. For a full
    # periodic/aperiodic decomposition use the `specparam` (FOOOF) package.
    fit = (freqs >= 2) & (freqs <= 40) & ~((freqs >= 7) & (freqs <= 14))
    slope, offset = np.polyfit(np.log10(freqs[fit]), np.log10(psd[:, fit].mean(0)), 1)
    feats["aperiodic_exponent"] = float(-slope)
    feats["aperiodic_offset"] = float(offset)
    return feats


def wpli_matrix(epochs: mne.Epochs, fmin: float, fmax: float) -> np.ndarray:
    """Weighted phase-lag index (Vinck et al., 2011), channels x channels."""
    data = epochs.get_data(copy=False)  # (epochs, channels, times)
    n_times = data.shape[-1]
    freqs = np.fft.rfftfreq(n_times, 1 / epochs.info["sfreq"])
    mask = (freqs >= fmin) & (freqs < fmax)
    taper = np.hanning(n_times)
    n_ch = data.shape[1]
    num = np.zeros((n_ch, n_ch, mask.sum()))
    den = np.zeros_like(num)
    for ep in data:  # accumulate per epoch to keep memory flat for dense nets
        X = np.fft.rfft(ep * taper, axis=-1)[:, mask]
        im = np.imag(X[:, None, :] * np.conj(X[None, :, :]))
        num += im
        den += np.abs(im)
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.abs(num) / den
    w = np.nan_to_num(w).mean(axis=-1)
    np.fill_diagonal(w, 0.0)
    return w


def connectivity_features(epochs: mne.Epochs, regions: dict[str, str]) -> dict[str, float]:
    feats: dict[str, float] = {}
    n_ch = len(epochs.ch_names)
    off_diag = ~np.eye(n_ch, dtype=bool)
    for band, (lo, hi) in BANDS.items():
        W = wpli_matrix(epochs, lo, hi)
        feats[f"wpli_{band}_global"] = float(W[off_diag].mean())
        strength = W.sum(axis=1) / (n_ch - 1)
        for region in REGIONS:
            idx = [i for i, ch in enumerate(epochs.ch_names) if regions.get(ch) == region]
            if idx:
                feats[f"wpli_{band}_{region}"] = float(strength[idx].mean())
    return feats


def extract_features(
    raw: mne.io.BaseRaw | str | Path,
    line_freq: float | None = 50.0,
    epoch_length: float = 2.0,
    reject_uv: float = 150.0,
) -> pd.Series:
    if not isinstance(raw, mne.io.BaseRaw):
        raw = load_raw(raw)
    raw = preprocess(raw, line_freq=line_freq)
    epochs = make_epochs(raw, epoch_length, reject_uv)
    regions = channel_regions(epochs.info)
    # qc_ columns are kept for reporting and dropped by data.load_cohort
    feats = {"qc_n_clean_epochs": float(len(epochs))}
    feats.update(spectral_features(epochs, regions))
    feats.update(connectivity_features(epochs, regions))
    return pd.Series(feats)
