#!/usr/bin/env python3
"""CLI entry point for extracting one or more HMM profiles by name (or accession) from an HMM
database such as Pfam-A.hmm, writing one <name>.hmm file per profile (ASCII HMMER3 format).
"""

import argparse
import sys
from pathlib import Path

from pyhmmer.plan7 import HMMFile


def _text(value: bytes | str | None) -> str:
    if value is None:
        return ""
    return value.decode() if isinstance(value, bytes) else value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("names", nargs="+",
                        help="Profile NAME(s) (e.g. DUF6365) or accession(s) (e.g. PF19892, version optional).")
    parser.add_argument("--db", type=Path, required=True, help="HMM database, e.g. data/pfam/Pfam-A.hmm.")
    parser.add_argument("--out-dir", type=Path, default=Path("."),
                        help="Directory to write one <name>.hmm file per profile into.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.db.exists():
        print(f"error: {args.db} not found", file=sys.stderr)
        return 1

    wanted = set(args.names)
    found: set[str] = set()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with HMMFile(args.db) as hmm_file:
        for hmm in hmm_file:
            name = _text(hmm.name)
            accession = _text(hmm.accession)
            keys = {name, accession, accession.split(".")[0]} & wanted
            if keys:
                out_path = args.out_dir / f"{name}.hmm"
                with out_path.open("wb") as fh:
                    hmm.write(fh, binary=False)
                print(f"{name} ({accession}) -> {out_path}")
                found |= keys
                if found == wanted:
                    break

    missing = sorted(wanted - found)
    if missing:
        print(f"WARNING: not found in {args.db}: {', '.join(missing)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
