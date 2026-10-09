"""Make the `gutbrain` package importable when scripts run from a checkout."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
