"""Loading curated product names from data/subcluster_groups.json, for RULE DESCRIPTION fields.

The json's own "id" field uses an underscore ("SCG_0031") where segment-dir/rule names don't
("SCG0031") -- normalised away on load so callers can look up by the plain SCG id.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

# antiSMASH's rule grammar tokenises DESCRIPTION text itself (rule_parser.py's Tokeniser):
# a character only continues a token if it's alphanumeric, '-', '_' or '.'; anything else must be
# one of the mapped structural characters ('(', ')', ',', whitespace) or the tokeniser raises
# RuleSyntaxError. Curated product names can contain apostrophes (e.g. "4'-O-carbamoyl-..."),
# which aren't in either set, so they must be stripped before use in a RULE block.
_UNSAFE_DESCRIPTION_CHARS = re.compile(r"[^\w\s().,-]")


def sanitize_description(name: str) -> str:
    """Strip characters antiSMASH's rule tokeniser can't handle from a DESCRIPTION value."""
    return _UNSAFE_DESCRIPTION_CHARS.sub("", name)


def load_product_names(path: Path) -> dict[str, str]:
    """SCG id (e.g. 'SCG0031') -> curated product name (e.g. '3,5-dihydroxyphenylglycine')."""
    data = json.loads(Path(path).read_text())
    return {entry["id"].replace("_", ""): entry["product"]["name"] for entry in data}
