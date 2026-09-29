#!/usr/bin/env python3
"""CLI entry point for Pfam domain mapping (docs/segment_grouping_spec.md, Phase 2).

Reads phase-1 outputs (segments.tsv / groups.tsv, from scripts/segment_groups.py) and maps Pfam
domains -- found by hmmsearch on the original full-length proteins -- onto them. --fasta-dir has
the same single-group-or-parent-directory semantics as scripts/segment_groups.py.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.pfam.params import PfamParams
from rule_suggester.pfam.pipeline import run_pfam

DATA_PFAM_DIR = Path(__file__).resolve().parent.parent / "data" / "pfam"
DEFAULT_PFAM_HMM = DATA_PFAM_DIR / "Pfam-A.hmm"
DEFAULT_PFAM_CLANS = DATA_PFAM_DIR / "Pfam-A.clans.tsv"


def parse_args() -> argparse.Namespace:
    defaults = PfamParams(hmm_path=DEFAULT_PFAM_HMM)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta-dir", type=Path, required=True,
                         help="A single group directory, or a parent directory of group directories.")
    parser.add_argument("--segments-dir", type=Path, required=True,
                         help="Phase-1 output root (the --out-dir given to segment_groups.py).")
    parser.add_argument("--pfam-hmm", type=Path, default=DEFAULT_PFAM_HMM)
    parser.add_argument("--pfam-clans", type=Path, default=DEFAULT_PFAM_CLANS,
                         help="Pfam-A.clans.tsv, for reporting each signature domain's clan.")
    parser.add_argument("--min-prevalence", type=float, default=defaults.min_prevalence)
    parser.add_argument("--straddle-threshold", type=float, default=defaults.straddle_threshold)
    parser.add_argument("--cpus", type=int, default=defaults.cpus)
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.pfam_hmm.exists():
        print(f"Pfam HMM file not found at {args.pfam_hmm}. Stopping.", file=sys.stderr)
        return 1
    if not args.pfam_clans.exists():
        print(f"Pfam clans file not found at {args.pfam_clans}. Stopping.", file=sys.stderr)
        return 1

    params = PfamParams(
        hmm_path=args.pfam_hmm,
        min_prevalence=args.min_prevalence,
        straddle_threshold=args.straddle_threshold,
        cpus=args.cpus,
    )

    try:
        results = run_pfam(args.fasta_dir, args.segments_dir, params, clans_path=args.pfam_clans)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    print(f"\nDone: {len(results)} group(s) processed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
