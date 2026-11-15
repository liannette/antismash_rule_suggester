#!/usr/bin/env python3
"""CLI entry point for extracting per-domain HMM profiles for every identifier referenced in a
combined rules file (scripts/combine_rules.py's output), from Pfam-A.hmm and the NCBI PGAP library.
PGAP models are written under their rule name (rule_suggester.pfam.libraries.pgap_rule_name), the
same name step 2 gave them, with the HMM NAME line set to it.

Writes one <domain>.hmm file per referenced profile (ASCII HMMER3 format), matching the
one-file-per-domain convention antiSMASH's subclusters module data/ directory already uses
(antismash/antismash/detection/subclusters/data/) -- ready to drop straight into it, alongside a
matching hmmdetails.txt entry (not generated here).
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.export.parse import extract_domain_ids
from rule_suggester.pfam.libraries import PgapLibrary, iter_hmm_file, load_pfam_names


def _text(value: bytes | str) -> str:
    return value.decode() if isinstance(value, bytes) else value

DATA = Path(__file__).resolve().parent.parent / "data"
DEFAULT_PFAM_HMM = DATA / "pfam" / "Pfam-A.hmm"
DEFAULT_PFAM_CLANS = DATA / "pfam" / "Pfam-A.clans.tsv"
DEFAULT_PGAP_HMM = DATA / "ncbi_pgap" / "hmm_PGAP.HMM.tgz"
DEFAULT_PGAP_METADATA = DATA / "ncbi_pgap" / "hmm_PGAP.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rules-file", type=Path, required=True,
                         help="A combined rules file, e.g. output/subcluster_rules.txt.")
    parser.add_argument("--pfam-hmm", type=Path, default=DEFAULT_PFAM_HMM)
    parser.add_argument("--pfam-clans", type=Path, default=DEFAULT_PFAM_CLANS,
                         help="Pfam-A.clans.tsv, for the Pfam names PGAP rule names must not collide with.")
    parser.add_argument("--pgap-hmm", type=Path, default=DEFAULT_PGAP_HMM)
    parser.add_argument("--pgap-metadata", type=Path, default=DEFAULT_PGAP_METADATA)
    parser.add_argument("--no-pgap", action="store_true", help="Look in Pfam only.")
    parser.add_argument("--out-dir", type=Path, required=True,
                         help="Directory to write one <domain>.hmm file per referenced profile into.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.rules_file.exists():
        print(f"error: {args.rules_file} not found", file=sys.stderr)
        return 1
    needed = [args.pfam_hmm] + ([] if args.no_pgap else [args.pfam_clans, args.pgap_hmm, args.pgap_metadata])
    for path in needed:
        if not path.exists():
            print(f"error: {path} not found", file=sys.stderr)
            return 1

    wanted = extract_domain_ids(args.rules_file.read_text())
    print(f"{len(wanted)} distinct domain(s) referenced in {args.rules_file}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    found: set[str] = set()
    libraries = [(args.pfam_hmm.name, lambda: iter_hmm_file(args.pfam_hmm))]
    if not args.no_pgap:
        pgap = PgapLibrary.load(args.pgap_hmm, args.pgap_metadata, load_pfam_names(args.pfam_clans))
        libraries.append((args.pgap_hmm.name, pgap.iter_hmms))
    for label, hmms in libraries:
        if found == wanted:
            break
        print(f"Scanning {label} ...")
        for hmm in hmms():
            name = _text(hmm.name)
            if name in wanted and name not in found:
                out_path = args.out_dir / f"{name}.hmm"
                with out_path.open("wb") as fh:
                    hmm.write(fh, binary=False)
                found.add(name)
                if found == wanted:
                    break

    missing = sorted(wanted - found)
    print(f"{len(found)}/{len(wanted)} written to {args.out_dir}")
    if missing:
        print(f"WARNING: {len(missing)} domain(s) not found in {' / '.join(l for l, _ in libraries)}:", file=sys.stderr)
        for name in missing:
            print(f"  {name}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
