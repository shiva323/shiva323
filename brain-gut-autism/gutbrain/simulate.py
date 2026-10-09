"""Synthetic ASD / TD cohort with a *known* gut -> brain -> behaviour pathway.

Purpose: run and debug the whole pipeline before real data arrive, and check
that each method recovers what was planted. Nothing here is a finding.

Planted structure (directions loosely follow the published literature):

* ASD children carry a latent dysbiosis factor that lowers Bifidobacterium,
  Prevotella, Coprococcus, Faecalibacterium and raises Desulfovibrio,
  Clostridium, Sutterella.
* Diet (fibre score) is lower in ASD and independently shapes the microbiome,
  i.e. diet is a confounder (Yap et al., Cell 2021).
* gut score = z(Bifidobacterium) - z(Desulfovibrio) drives frontal alpha wPLI
  and the aperiodic exponent of the EEG.
* SRS (social responsiveness) depends on diagnosis and on those two EEG
  features -> mediation Bifidobacterium -> frontal alpha wPLI -> SRS.
* Fecal metabolites (butyrate, kynurenate, p-cresol, ...) follow their
  producer taxa.
* ~15 % of children have no EEG and ~25 % no metabolomics (missing blocks).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .data import Cohort
from .eeg import BANDS, REGIONS

KEY_TAXA = {
    # genus: (dysbiosis loading, diet loading, baseline log-abundance)
    "Bifidobacterium": (-0.9, 0.3, 2.5),
    "Prevotella": (-0.6, 0.7, 2.0),
    "Coprococcus": (-0.5, 0.2, 1.5),
    "Faecalibacterium": (-0.4, 0.5, 3.0),
    "Roseburia": (-0.1, 0.5, 2.0),
    "Desulfovibrio": (0.8, 0.0, 0.5),
    "Clostridium": (0.6, -0.1, 1.5),
    "Sutterella": (0.5, 0.0, 0.8),
    "Bacteroides": (0.2, -0.3, 3.5),
    "Akkermansia": (0.0, 0.2, 1.0),
    "Lactobacillus": (0.0, 0.1, 0.5),
    "Ruminococcus": (0.1, 0.1, 2.0),
    "Blautia": (0.0, 0.0, 2.5),
    "Parabacteroides": (0.2, 0.0, 1.5),
    "Veillonella": (0.0, 0.0, 0.5),
}

METABOLITES = {
    # metabolite: {producer taxon: loading}, dysbiosis loading
    "butyrate": ({"Faecalibacterium": 0.5, "Roseburia": 0.4, "Coprococcus": 0.3}, 0.0),
    "propionate": ({"Bacteroides": 0.4, "Akkermansia": 0.3}, 0.0),
    "kynurenate": ({"Bifidobacterium": 0.3}, -0.5),
    "indole_3_propionate": ({"Clostridium": -0.2, "Coprococcus": 0.3}, -0.3),
    "p_cresol_sulfate": ({"Clostridium": 0.6}, 0.2),
    "hydrogen_sulfide_proxy": ({"Desulfovibrio": 0.7}, 0.0),
    "serotonin": ({"Clostridium": 0.2}, 0.3),
    "GABA": ({"Bifidobacterium": 0.4, "Lactobacillus": 0.3}, 0.0),
    "tryptophan": ({}, 0.0),
    "4_ethylphenyl_sulfate": ({"Clostridium": 0.3}, 0.3),
}


def _z(x: np.ndarray) -> np.ndarray:
    return (x - x.mean()) / x.std()


def simulate_cohort(
    n_asd: int = 120,
    n_td: int = 100,
    n_taxa: int = 120,
    n_noise_metabolites: int = 20,
    frac_missing_brain: float = 0.15,
    frac_missing_metab: float = 0.25,
    seed: int = 7,
) -> Cohort:
    rng = np.random.default_rng(seed)
    n = n_asd + n_td
    ids = [f"sub-{i:04d}" for i in range(1, n + 1)]
    asd = np.r_[np.ones(n_asd), np.zeros(n_td)]

    # ---- demographics, diet, site -------------------------------------------
    age = rng.uniform(3, 12, n)
    sex = np.where(rng.random(n) < np.where(asd == 1, 0.8, 0.55), "M", "F")
    fiber = rng.normal(-0.6 * asd, 1.0)  # selective eating in ASD
    site = rng.choice(["site_A", "site_B"], n)
    dysbiosis = 1.3 * asd + rng.normal(0, 1, n)

    # ---- microbiome counts ----------------------------------------------------
    names = list(KEY_TAXA) + [f"Genus_{i:03d}" for i in range(n_taxa - len(KEY_TAXA))]
    base = np.r_[[v[2] for v in KEY_TAXA.values()], rng.normal(-0.5, 1.5, n_taxa - len(KEY_TAXA))]
    load_dys = np.r_[[v[0] for v in KEY_TAXA.values()], np.zeros(n_taxa - len(KEY_TAXA))]
    load_diet = np.r_[[v[1] for v in KEY_TAXA.values()], rng.normal(0, 0.15, n_taxa - len(KEY_TAXA))]
    batch = np.zeros(n_taxa)
    batch[rng.choice(n_taxa, 20, replace=False)] = rng.normal(0, 0.6, 20)
    logab = (
        base
        + np.outer(dysbiosis, load_dys)
        + np.outer(fiber, load_diet)
        + np.outer(site == "site_B", batch)
        + rng.normal(0, 0.8, (n, n_taxa))
    )
    props = np.exp(logab)
    props /= props.sum(axis=1, keepdims=True)
    depth = rng.lognormal(np.log(30000), 0.4, n).astype(int)
    counts = np.array([rng.multinomial(d, p) for d, p in zip(depth, props)])
    microbiome = pd.DataFrame(counts, index=ids, columns=names)

    clr = np.log(counts + 0.5)
    clr -= clr.mean(axis=1, keepdims=True)
    col = {name: clr[:, i] for i, name in enumerate(names)}
    gut_score = _z(_z(col["Bifidobacterium"]) - _z(col["Desulfovibrio"]))

    # ---- EEG-like brain features ---------------------------------------------
    brain = {}
    band_base = np.array([1.2, 0.6, 0.9, 0.3, -0.8])  # delta..gamma log power
    for region in REGIONS:
        logp = band_base + rng.normal(0, 0.25, (n, 5))
        logp[:, 0] -= 0.06 * (age - 7)  # delta falls with age
        logp[:, 2] += 0.05 * (age - 7) + (0.25 if region == "occipital" else 0)
        if region == "temporal":
            logp[:, 4] += 0.25 * asd  # higher temporal gamma in ASD
        rel = np.exp(logp) / np.exp(logp).sum(axis=1, keepdims=True)
        for b, band in enumerate(BANDS):
            brain[f"relpow_{band}_{region}"] = rel[:, b]
    brain["theta_beta_ratio_frontal"] = np.exp(rng.normal(0.9 + 0.1 * asd, 0.25))
    brain["alpha_peak_freq"] = 8.0 + 0.2 * (age - 3) + rng.normal(0, 0.6, n)
    brain["aperiodic_exponent"] = (
        1.6 - 0.03 * (age - 7) + 0.12 * (0.6 * gut_score + 0.8 * rng.normal(0, 1, n))
    )
    brain["aperiodic_offset"] = rng.normal(1.0, 0.2, n)
    for band in BANDS:
        glob = rng.normal(0, 1, n)
        brain[f"wpli_{band}_global"] = 0.15 + 0.03 * glob
        for region in REGIONS:
            brain[f"wpli_{band}_{region}"] = 0.15 + 0.03 * (0.6 * glob + 0.8 * rng.normal(0, 1, n))
    brain["wpli_alpha_frontal"] = 0.15 + 0.03 * (0.55 * gut_score + 0.83 * rng.normal(0, 1, n))
    brain = pd.DataFrame(brain, index=ids)

    # ---- metabolites (log2 intensities) --------------------------------------
    metab = {}
    for met, (producers, dys) in METABOLITES.items():
        v = dys * dysbiosis + rng.normal(0, 0.8, n)
        for taxon, w in producers.items():
            v += w * _z(col[taxon])
        metab[met] = 15 + v
    for i in range(n_noise_metabolites):
        metab[f"metabolite_{i:02d}"] = 15 + rng.normal(0, 1, n)
    metab = pd.DataFrame(metab, index=ids)

    # ---- clinical scores ------------------------------------------------------
    srs = (
        50
        + 25 * asd
        - 4.0 * _z(brain["wpli_alpha_frontal"].to_numpy())
        - 2.5 * _z(brain["aperiodic_exponent"].to_numpy())
        + rng.normal(0, 7, n)
    )
    gsrs = np.clip(8 + 4 * asd + 2.5 * _z(col["Desulfovibrio"]) + rng.normal(0, 3, n), 0, None)
    ados = np.where(asd == 1, np.clip(np.round(6 + 1.0 * rng.normal(0, 1.5, n) - 0.8 * gut_score), 1, 10), np.nan)

    subjects = pd.DataFrame(
        {
            "group": np.where(asd == 1, "ASD", "TD"),
            "age": age.round(1),
            "sex": sex,
            "site": site,
            "fiber_score": fiber.round(2),
            "SRS_total": np.clip(srs, 30, 90).round(0),
            "ADOS_CSS": ados,
            "GSRS_total": gsrs.round(1),
        },
        index=ids,
    )

    keep_brain = rng.random(n) >= frac_missing_brain
    keep_metab = rng.random(n) >= frac_missing_metab
    return Cohort(
        subjects=subjects,
        microbiome=microbiome,
        brain=brain.loc[keep_brain],
        metabolites=metab.loc[keep_metab],
    )
