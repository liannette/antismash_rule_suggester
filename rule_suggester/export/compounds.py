"""Loading curated compound metadata from data/building_blocks.json (subcluster_collection's
build_jsons.py output), for building compound_details.txt (antiSMASH subclusters module format --
see antismash/antismash/detection/subclusters/compounds.py).
"""

from __future__ import annotations

import json
from pathlib import Path


def load_products(path: Path) -> dict[str, dict]:
    """bb_id (e.g. 'BB0003') -> {"name", "smiles", "classification"} for that building block.
    classification is [bb_chemical_class] (one controlled-vocabulary value, or empty)."""
    data = json.loads(Path(path).read_text())
    return {
        entry["id"]: {
            "name": entry.get("name"),
            "smiles": entry.get("smiles"),
            "classification": [entry["chemical_class"]] if entry.get("chemical_class") else [],
        }
        for entry in data
    }
