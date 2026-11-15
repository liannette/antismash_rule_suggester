"""Loading curated building-block names from data/building_blocks.json (subcluster_collection's
build_jsons.py output), for RULE DESCRIPTION fields. Keyed by bb_id, which is also the prefix of
each segment/rule directory name (BB0003_L-dihydrostreptose -> RULE BB0003).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

# antiSMASH's rule grammar tokenises DESCRIPTION text itself (rule_parser.py's Tokeniser):
# a character only continues a token if it's alphanumeric, '-', '_' or '.'; anything else must be
# one of the mapped structural characters ('(', ')', ',', whitespace) or the tokeniser raises
# RuleSyntaxError. Curated names can contain apostrophes (e.g. "4'-O-carbamoyl-..."),
# which aren't in either set, so they must be stripped before use in a RULE block.
_UNSAFE_DESCRIPTION_CHARS = re.compile(r"[^\w\s().,-]")


def sanitize_description(name: str) -> str:
    """Strip characters antiSMASH's rule tokeniser can't handle from a DESCRIPTION value."""
    return _UNSAFE_DESCRIPTION_CHARS.sub("", name)


def load_product_names(path: Path) -> dict[str, str]:
    """bb_id (e.g. 'BB0003') -> curated building-block name (e.g. 'L-dihydrostreptose')."""
    data = json.loads(Path(path).read_text())
    return {entry["id"]: entry["name"] for entry in data if entry.get("name")}
