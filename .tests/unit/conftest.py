"""Make the workflow scripts importable from the unit tests."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "workflow" / "scripts"))
