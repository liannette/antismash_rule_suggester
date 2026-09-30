#!/usr/bin/env python3
"""CLI entry point for combining a segment group's report.md (Phase 1), pfam_report.md
(Phase 2) and rule.txt (Phase 3) into one combined_report.md.

--segments-dir accepts the same single-group-or-parent-directory semantics as the other
scripts: point it at one output/segments/<group> directory to process just that group, or at
output/segments itself to process every group directory found inside it.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.report import write_combined_report


def _is_group_dir(path: Path) -> bool:
    return (path / "report.md").exists()


def discover_group_dirs(segments_dir: Path) -> list[Path]:
    if not segments_dir.is_dir():
        raise ValueError(f"--segments-dir {segments_dir} is not a directory")
    if _is_group_dir(segments_dir):
        return [segments_dir]
    group_dirs = sorted(p for p in segments_dir.iterdir() if p.is_dir() and _is_group_dir(p))
    if not group_dirs:
        raise ValueError(f"{segments_dir} contains no group directories with a report.md")
    return group_dirs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segments-dir", type=Path, required=True,
                         help="A single group directory, or a parent directory of group directories.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        group_dirs = discover_group_dirs(args.segments_dir)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    for group_dir in group_dirs:
        out_path = write_combined_report(group_dir)
        print(f"  {group_dir.name}: {out_path}")

    print(f"Done: {len(group_dirs)} group(s) processed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
