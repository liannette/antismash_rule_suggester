"""Loading curated compound metadata from data/subcluster_groups.json, for building
compound_details.txt (antiSMASH subclusters module format -- see
antismash/antismash/detection/subclusters/compounds.py).
"""

from __future__ import annotations

import json
from pathlib import Path


def load_products(path: Path) -> dict[str, dict]:
    """SCG id (e.g. 'SCG0031') -> that entry's curated "product" dict (name/smiles/classification/...)."""
    data = json.loads(Path(path).read_text())
    return {entry["id"].replace("_", ""): entry["product"] for entry in data}
