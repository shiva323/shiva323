"""Write a synthetic ASD/TD gut–brain cohort in the pipeline's CSV layout.

    python scripts/make_synthetic.py --out data/synthetic
"""
import argparse

import _path  # noqa: F401

from gutbrain.data import save_cohort
from gutbrain.simulate import simulate_cohort

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--out", default="data/synthetic")
parser.add_argument("--n-asd", type=int, default=120)
parser.add_argument("--n-td", type=int, default=100)
parser.add_argument("--seed", type=int, default=7)
args = parser.parse_args()

cohort = simulate_cohort(n_asd=args.n_asd, n_td=args.n_td, seed=args.seed)
save_cohort(cohort, args.out)
print(cohort.summary())
print(f"written to {args.out}/")
