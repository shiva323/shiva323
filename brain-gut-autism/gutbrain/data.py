"""Load and align the per-subject tables of a gut–brain cohort.

Expected layout of a cohort directory (every table is keyed by ``subject_id``)::

    subjects.csv            required  group (ASD/TD), age, sex + any clinical scores
    microbiome_counts.csv   required  raw read counts, one column per taxon (e.g. genus)
    brain_features.csv      optional  EEG / sMRI / fMRI derived features
    metabolites.csv         optional  fecal or plasma metabolite intensities

Subjects may be missing from the optional tables (a child without a usable
EEG, for instance); the late-fusion model handles that case.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

GROUP_POSITIVE = "ASD"
GROUP_NEGATIVE = "TD"


@dataclass
class Cohort:
    subjects: pd.DataFrame
    microbiome: pd.DataFrame
    brain: pd.DataFrame | None = None
    metabolites: pd.DataFrame | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def y(self) -> pd.Series:
        """Binary diagnosis label (1 = ASD, 0 = TD) indexed by subject."""
        return (self.subjects["group"] == GROUP_POSITIVE).astype(int)

    def covariates(self, names: list[str]) -> pd.DataFrame:
        """Numeric covariate matrix; categorical columns are one-hot encoded."""
        missing = [c for c in names if c not in self.subjects.columns]
        if missing:
            raise KeyError(f"covariates not found in subjects.csv: {missing}")
        cov = pd.get_dummies(self.subjects[names], drop_first=True, dtype=float)
        if cov.isna().any().any():
            # Median imputation keeps every subject; report it so it is not silent.
            n = int(cov.isna().sum().sum())
            self.notes.append(f"imputed {n} missing covariate values with the median")
            cov = cov.fillna(cov.median())
        return cov.astype(float)

    def modalities(self) -> dict[str, pd.DataFrame]:
        out = {"microbiome": self.microbiome}
        if self.brain is not None:
            out["brain"] = self.brain
        if self.metabolites is not None:
            out["metabolites"] = self.metabolites
        return out

    def summary(self) -> str:
        counts = self.subjects["group"].value_counts().to_dict()
        lines = [f"subjects: {len(self.subjects)} {counts}"]
        for name, df in self.modalities().items():
            lines.append(f"{name}: {df.shape[0]} subjects x {df.shape[1]} features")
        return "\n".join(lines)


def _read(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "subject_id" not in df.columns:
        raise ValueError(f"{path.name}: missing 'subject_id' column")
    if df["subject_id"].duplicated().any():
        dup = df.loc[df["subject_id"].duplicated(), "subject_id"].tolist()[:5]
        raise ValueError(f"{path.name}: duplicated subject_id values, e.g. {dup}")
    return df.set_index("subject_id")


def load_cohort(directory: str | Path) -> Cohort:
    directory = Path(directory)
    subjects = _read(directory / "subjects.csv")
    if "group" not in subjects.columns:
        raise ValueError("subjects.csv: missing 'group' column (values ASD / TD)")
    bad = set(subjects["group"].unique()) - {GROUP_POSITIVE, GROUP_NEGATIVE}
    if bad:
        raise ValueError(f"subjects.csv: unexpected group labels {bad}; use ASD / TD")

    microbiome = _read(directory / "microbiome_counts.csv")
    if (microbiome.values < 0).any():
        raise ValueError("microbiome_counts.csv must hold non-negative read counts")

    def optional(name: str) -> pd.DataFrame | None:
        path = directory / name
        return _read(path).astype(float) if path.exists() else None

    cohort = Cohort(
        subjects=subjects,
        microbiome=microbiome.astype(float),
        brain=optional("brain_features.csv"),
        metabolites=optional("metabolites.csv"),
    )

    # Every analysis is anchored on subjects that have microbiome data and a label.
    keep = subjects.index.intersection(microbiome.index)
    dropped = len(subjects) - len(keep)
    if dropped:
        cohort.notes.append(f"dropped {dropped} subjects without microbiome data")
    cohort.subjects = subjects.loc[keep]
    cohort.microbiome = cohort.microbiome.loc[keep]
    for attr in ("brain", "metabolites"):
        df = getattr(cohort, attr)
        if df is not None:
            df = df.loc[df.index.intersection(keep)]
            df = df.loc[:, ~df.columns.str.startswith("qc_")]  # quality-control only
            df = df.loc[:, df.notna().mean() > 0.5]  # drop mostly-empty features
            if df.isna().any().any():
                cohort.notes.append(f"{attr}: median-imputed sparse missing values")
                df = df.fillna(df.median())
            setattr(cohort, attr, df)
    return cohort


def save_cohort(cohort: Cohort, directory: str | Path) -> None:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    cohort.subjects.to_csv(directory / "subjects.csv", index_label="subject_id")
    cohort.microbiome.round().astype(np.int64).to_csv(
        directory / "microbiome_counts.csv", index_label="subject_id"
    )
    if cohort.brain is not None:
        cohort.brain.to_csv(directory / "brain_features.csv", index_label="subject_id")
    if cohort.metabolites is not None:
        cohort.metabolites.to_csv(directory / "metabolites.csv", index_label="subject_id")
