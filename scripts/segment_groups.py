#!/usr/bin/env python3
"""CLI entry point for the segment-grouping pipeline (docs/segment_grouping_spec.md).

--fasta-dir accepts either a single product-group directory or a parent directory containing
several of them; which case applies is detected from its contents.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.segments.io_fasta import discover_group_dirs, is_group_dir
from rule_suggester.segments.params import Params
from rule_suggester.segments.pipeline import run_group

DEFAULT_DIAMOND_BIN = Path("/Users/alien/miniconda3/envs/rule_suggester/bin/diamond")


def parse_args() -> argparse.Namespace:
    defaults = Params()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta-dir", type=Path, required=True,
                         help="A single group directory, or a parent directory of group directories.")
    parser.add_argument("--out-dir", type=Path, required=True,
                         help="Output root; each group's outputs go in <out-dir>/<group-name>/.")
    parser.add_argument("--diamond-bin", type=Path, default=DEFAULT_DIAMOND_BIN)
    parser.add_argument("--evalue", type=float, default=defaults.evalue)
    parser.add_argument("--min-hsp-len", type=int, default=defaults.min_hsp_len)
    parser.add_argument("--tol", type=int, default=defaults.tol)
    parser.add_argument("--min-seg", type=int, default=defaults.min_seg)
    parser.add_argument("--min-recip-overlap", type=float, default=defaults.min_recip_overlap)
    parser.add_argument("--min-w", type=float, default=defaults.min_w)
    parser.add_argument("--min-link", type=float, default=defaults.min_link)
    parser.add_argument("--max-iter", type=int, default=defaults.max_iter)
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.diamond_bin.exists():
        print(f"DIAMOND not found at {args.diamond_bin}. Stopping.", file=sys.stderr)
        return 1

    params = Params(
        evalue=args.evalue,
        min_hsp_len=args.min_hsp_len,
        tol=args.tol,
        min_seg=args.min_seg,
        min_recip_overlap=args.min_recip_overlap,
        min_w=args.min_w,
        min_link=args.min_link,
        max_iter=args.max_iter,
    )

    try:
        group_dirs = discover_group_dirs(args.fasta_dir)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    single_group = is_group_dir(args.fasta_dir)

    failed = []
    for group_dir in group_dirs:
        out_dir = args.out_dir if single_group else args.out_dir / group_dir.name
        print(f"Processing {group_dir.name} ...")
        result = run_group(group_dir, out_dir, params, args.diamond_bin)
        status = "OK" if not result.invariant_violations else f"{len(result.invariant_violations)} VIOLATION(S)"
        print(
            f"  {result.n_scs} subclusters, {result.n_proteins} proteins, {result.n_segments} segments, "
            f"{result.n_groups} groups -> {out_dir} [{status}]"
        )
        counts = ", ".join(f"{gid}: {n}" for gid, n in result.group_sc_counts)
        print(f"    subcluster instances per group -- {counts}")
        if result.invariant_violations:
            failed.append(result)
            for v in result.invariant_violations:
                print(f"    - {v}", file=sys.stderr)

    if failed:
        print(f"\n{len(failed)}/{len(group_dirs)} group(s) had invariant violations.", file=sys.stderr)
        return 1

    print(f"\nDone: {len(group_dirs)} group(s) processed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
