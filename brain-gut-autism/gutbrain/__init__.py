"""gutbrain: multimodal gut–brain axis analysis toolkit for autism (ASD) cohorts.

Modules
-------
data         load / validate / align the per-subject CSV tables
simulate     synthetic cohort with a known gut -> brain -> behaviour pathway
microbiome   filtering, CLR transform, alpha / beta diversity, PERMANOVA
eeg          raw EEG -> spectral, aperiodic and wPLI connectivity features
fmri         ROI time series -> functional connectivity features
stats        covariate-adjusted group tests (FDR) and mediation analysis
integration  sparse CCA between the gut and brain blocks
models       nested-CV classifiers: single modality, early and late fusion
plots        figures used by the report
"""

__version__ = "0.1.0"
