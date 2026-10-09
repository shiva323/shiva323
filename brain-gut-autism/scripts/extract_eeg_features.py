"""Extract resting-state EEG features for many subjects -> brain_features.csv

    python scripts/extract_eeg_features.py --files "eeg/sub-*_task-rest_eeg.set" \
        --id-regex "(sub-[A-Za-z0-9]+)" --line-freq 50 --out data/mycohort/brain_features.csv

Each file is one subject; the subject_id is taken from the file name with
--id-regex and must match subject_id in subjects.csv.
"""
import argparse
import glob
import re

import _path  # noqa: F401
import pandas as pd

from gutbrain.eeg import extract_features

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--files", required=True, help="glob pattern of EEG recordings")
parser.add_argument("--id-regex", default=r"(sub-[A-Za-z0-9]+)")
parser.add_argument("--line-freq", type=float, default=50.0, help="50 (Iran/Europe) or 60 (US)")
parser.add_argument("--epoch-length", type=float, default=2.0)
parser.add_argument("--reject-uv", type=float, default=150.0)
parser.add_argument("--out", required=True)
args = parser.parse_args()

rows = {}
for path in sorted(glob.glob(args.files)):
    match = re.search(args.id_regex, path)
    if not match:
        print(f"skip (no subject id): {path}")
        continue
    try:
        rows[match.group(1)] = extract_features(path, args.line_freq, args.epoch_length, args.reject_uv)
        print(f"ok   {match.group(1)}  epochs={int(rows[match.group(1)]['qc_n_clean_epochs'])}")
    except Exception as exc:  # keep going; one bad recording should not stop the batch
        print(f"FAIL {path}: {exc}")

if not rows:
    raise SystemExit("no recordings processed")
pd.DataFrame(rows).T.to_csv(args.out, index_label="subject_id")
print(f"{len(rows)} subjects -> {args.out}")
