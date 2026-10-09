"""Run the full gut–brain analysis on a cohort directory and write a report.

    python scripts/run_pipeline.py --data data/synthetic --out results/synthetic \
        --covariates age sex site fiber_score \
        --mediation Bifidobacterium:wpli_alpha_frontal:SRS_total

Cohort directory layout: see gutbrain/data.py (subjects.csv,
microbiome_counts.csv, optional brain_features.csv and metabolites.csv).
"""
import argparse

import _path  # noqa: F401

from gutbrain.data import load_cohort
from gutbrain.pipeline import PipelineConfig, run_pipeline

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--data", required=True, help="cohort directory with the CSV tables")
parser.add_argument("--out", required=True, help="output directory for report and tables")
parser.add_argument("--covariates", nargs="*", default=["age", "sex"],
                    help="columns of subjects.csv to adjust for (diet, site, BMI, meds ...)")
parser.add_argument("--mediation", nargs="*", default=[],
                    help="a-priori hypotheses written as X:M:Y, e.g. taxon:brain_feature:score")
parser.add_argument("--severity", default=None, help="severity score for the within-ASD regression")
parser.add_argument("--min-prevalence", type=float, default=0.1)
parser.add_argument("--n-perm", type=int, default=500, help="permutations for the sCCA test")
parser.add_argument("--cv-repeats", type=int, default=3)
parser.add_argument("--seed", type=int, default=0)
args = parser.parse_args()

mediation = []
for spec in args.mediation:
    parts = spec.split(":")
    if len(parts) != 3:
        parser.error(f"--mediation expects X:M:Y, got '{spec}'")
    mediation.append(tuple(parts))

cohort = load_cohort(args.data)
print(cohort.summary())
cfg = PipelineConfig(
    covariates=args.covariates,
    mediation=mediation,
    severity=args.severity,
    min_prevalence=args.min_prevalence,
    n_perm=args.n_perm,
    cv_repeats=args.cv_repeats,
    seed=args.seed,
)
summary = run_pipeline(cohort, args.out, cfg)
print(f"report: {args.out}/report.md")
for model, row in summary["classification"].items():
    print(f"  {model:13s} AUC = {row['AUC_mean']:.3f} ± {row['AUC_sd']:.3f}")
