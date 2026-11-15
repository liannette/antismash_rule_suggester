"""Extracting Pfam profile identifiers referenced in a combined rules file's CONDITIONS/EXTENDERS
lines -- for cross-referencing against a module's registered signatures, or for extracting the
matching HMM profiles out of a Pfam database.
"""

from __future__ import annotations

import re

_KEYWORDS = {
    "and", "or", "not", "cds", "repeatable", "minimum",
    "RULE", "CATEGORY", "DESCRIPTION", "CUTOFF", "NEIGHBOURHOOD", "CONDITIONS", "EXTENDERS",
    "Subcluster",
}
_IDENTIFIER_RE = re.compile(r"(?=[A-Za-z0-9_-]*[A-Za-z])[A-Za-z0-9_-]+")
_CONDITIONS_RE = re.compile(r"^\s*CONDITIONS(.*?)(?=^\s*EXTENDERS\b|\Z)", re.DOTALL | re.MULTILINE)
_EXTENDERS_RE = re.compile(r"^\s*EXTENDERS(.*?)\Z", re.DOTALL | re.MULTILINE)
_RULE_ID_RE = re.compile(r"^RULE\s+(\S+)", re.MULTILINE)
_COMMENT_RE = re.compile(r"#.*")


def strip_comments(rules_text: str) -> str:
    """Drop '#' comments (to end of line), e.g. the '# MANUAL' / '# WARNING' notes step 3 writes,
    so their words aren't read as identifiers or rule ids."""
    return _COMMENT_RE.sub("", rules_text)


def extract_rule_ids(rules_text: str) -> list[str]:
    """Every RULE id in a combined rules file, in file order -- e.g. to know which rows a
    matching compound_details.txt needs (antiSMASH KeyErrors on load if any RULE is missing one).
    """
    return _RULE_ID_RE.findall(strip_comments(rules_text))


def extract_domain_ids(rules_text: str) -> set[str]:
    """Every Pfam profile identifier referenced in any RULE block's CONDITIONS/EXTENDERS, across
    the whole file. Only those two fields are scanned -- not RULE/DESCRIPTION -- so a product
    name or rule id that happens to look like an identifier doesn't get pulled in by mistake.
    """
    ids: set[str] = set()
    for block in re.split(r"(?=^RULE )", strip_comments(rules_text), flags=re.MULTILINE):
        if not block.lstrip().startswith("RULE "):
            continue
        for pattern in (_CONDITIONS_RE, _EXTENDERS_RE):
            match = pattern.search(block)
            if match:
                ids.update(_IDENTIFIER_RE.findall(match.group(1)))
    return ids - _KEYWORDS
