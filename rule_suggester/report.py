"""Combining a segment group's per-phase reports (segment grouping, Pfam domains, rule
conditions) into one file, so reviewing a group doesn't mean opening report.md, pfam_report.md
and rule.txt separately.

Reads the already-written outputs back in rather than re-deriving anything -- each phase's own
script remains the source of truth for its content; this only stitches it together.
"""

from __future__ import annotations

import re
from pathlib import Path

_HEADING_RE = re.compile(r"^(#+)( .*)?$")
_SECTION_HEADING_RE = re.compile(r"^#+\s+(.*)$")
_WARNING_HEADINGS = ("Warnings", "Flags")  # Phase 1 calls its section "Flags"; 2 and 3 call it "Warnings"


def _demote_headings(text: str, by: int = 1) -> str:
    """Shift every Markdown heading down `by` levels, so a sub-report's own top-level title
    nests under this file's title instead of competing with it.
    """
    lines = []
    for line in text.splitlines():
        match = _HEADING_RE.match(line)
        lines.append(f"{'#' * by}{line}" if match else line)
    return "\n".join(lines)


def _extract_warnings(text: str) -> list[str]:
    """Pull the bullet/bold lines out of a sub-report's own Warnings/Flags section, so they can
    be rolled up to the top of the combined report -- unusual or dropped data three files deep
    is easy to miss otherwise.
    """
    lines = text.splitlines()
    collected: list[str] = []
    i = 0
    while i < len(lines):
        heading = _SECTION_HEADING_RE.match(lines[i])
        if heading and heading.group(1).strip() in _WARNING_HEADINGS:
            i += 1
            while i < len(lines) and not _SECTION_HEADING_RE.match(lines[i]):
                content = lines[i].strip()
                if content and content != "None.":
                    collected.append(content.lstrip("- "))
                i += 1
        else:
            i += 1
    return collected


def write_combined_report(group_dir: Path) -> Path:
    group_dir = Path(group_dir)
    group_name = group_dir.name

    steps = [
        ("Step 1: Segment grouping", group_dir / "report.md", "scripts/segment_groups.py"),
        ("Step 2: Pfam domains", group_dir / "pfam_report.md", "scripts/pfam_domains.py"),
        ("Step 3: Rule conditions", group_dir / "rule_conditions.md", "scripts/build_rule.py"),
    ]
    texts = [(label, path.read_text().strip() if path.exists() else None, path, script) for label, path, script in steps]

    lines = [f"# Subcluster report: {group_name}", ""]

    lines.append("## Warnings")
    lines.append("")
    warnings = [(label, w) for label, text, _, _ in texts if text for w in _extract_warnings(text)]
    if warnings:
        lines.extend(f"- [{label}] {w}" for label, w in warnings)
    else:
        lines.append("None.")
    lines.append("")

    for label, text, path, script in texts:
        lines.append(f"## {label}")
        lines.append("")
        if text is not None:
            lines.append(_demote_headings(text))
        else:
            lines.append(f"*{path.name} not found -- run {script} first.*")
        lines.append("")

    out_path = group_dir / "combined_report.md"
    out_path.write_text("\n".join(lines))
    return out_path
