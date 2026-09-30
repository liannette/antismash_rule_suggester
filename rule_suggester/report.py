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


def _demote_headings(text: str, by: int = 1) -> str:
    """Shift every Markdown heading down `by` levels, so a sub-report's own top-level title
    nests under this file's title instead of competing with it.
    """
    lines = []
    for line in text.splitlines():
        match = _HEADING_RE.match(line)
        lines.append(f"{'#' * by}{line}" if match else line)
    return "\n".join(lines)


def write_combined_report(group_dir: Path) -> Path:
    group_dir = Path(group_dir)
    group_name = group_dir.name
    lines = [f"# Subcluster report: {group_name}", ""]

    lines.append("## Step 1: Segment grouping")
    lines.append("")
    segment_report = group_dir / "report.md"
    if segment_report.exists():
        lines.append(_demote_headings(segment_report.read_text().strip()))
    else:
        lines.append(f"*{segment_report.name} not found -- run scripts/segment_groups.py first.*")
    lines.append("")

    lines.append("## Step 2: Pfam domains")
    lines.append("")
    pfam_report = group_dir / "pfam_report.md"
    if pfam_report.exists():
        lines.append(_demote_headings(pfam_report.read_text().strip()))
    else:
        lines.append(f"*{pfam_report.name} not found -- run scripts/pfam_domains.py first.*")
    lines.append("")

    lines.append("## Step 3: Rule conditions")
    lines.append("")
    rule_conditions = group_dir / "rule_conditions.md"
    if rule_conditions.exists():
        lines.append(_demote_headings(rule_conditions.read_text().strip()))
    else:
        lines.append(f"*{rule_conditions.name} not found -- run scripts/build_rule.py first.*")
    lines.append("")

    out_path = group_dir / "combined_report.md"
    out_path.write_text("\n".join(lines))
    return out_path
