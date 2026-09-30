#!/usr/bin/env python3
"""CLI entry point for extracting per-domain HMM profiles for every Pfam identifier referenced
in a combined rules file (scripts/combine_rules.py's output), from the full Pfam-A.hmm database.

Writes one <domain>.hmm file per referenced profile (ASCII HMMER3 format), matching the
one-file-per-domain convention antiSMASH's subclusters module data/ directory already uses
(antismash/antismash/detection/subclusters/data/) -- ready to drop straight into it, alongside a
matching hmmdetails.txt entry (not generated here).
"""

import argparse
import sys
from pathlib import Path

from pyhmmer.plan7 import HMMFile

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.export.parse import extract_domain_ids


def _text(value: bytes | str) -> str:
    return value.decode() if isinstance(value, bytes) else value

DEFAULT_PFAM_HMM = Path(__file__).resolve().parent.parent / "data" / "pfam" / "Pfam-A.hmm"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rules-file", type=Path, required=True,
                         help="A combined rules file, e.g. output/subcluster_rules.txt.")
    parser.add_argument("--pfam-hmm", type=Path, default=DEFAULT_PFAM_HMM)
    parser.add_argument("--out-dir", type=Path, required=True,
                         help="Directory to write one <domain>.hmm file per referenced profile into.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.rules_file.exists():
        print(f"error: {args.rules_file} not found", file=sys.stderr)
        return 1
    if not args.pfam_hmm.exists():
        print(f"error: {args.pfam_hmm} not found", file=sys.stderr)
        return 1

    wanted = extract_domain_ids(args.rules_file.read_text())
    print(f"{len(wanted)} distinct domain(s) referenced in {args.rules_file}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    found: set[str] = set()
    print(f"Scanning {args.pfam_hmm.name} ...")
    with HMMFile(args.pfam_hmm) as hmm_file:
        for hmm in hmm_file:
            name = _text(hmm.name)
            if name in wanted:
                out_path = args.out_dir / f"{name}.hmm"
                with out_path.open("wb") as fh:
                    hmm.write(fh, binary=False)
                found.add(name)
                if found == wanted:
                    break

    missing = sorted(wanted - found)
    print(f"{len(found)}/{len(wanted)} written to {args.out_dir}")
    if missing:
        print(f"WARNING: {len(missing)} domain(s) not found in {args.pfam_hmm}:", file=sys.stderr)
        for name in missing:
            print(f"  {name}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
