#!/usr/bin/env python3
"""CLI entry point for extracting per-domain HMM profiles for every identifier referenced in a
combined rules file (scripts/combine_rules.py's output), from Pfam-A.hmm and the NCBI PGAP library.
Curated smCOG profiles (data/smcogs/smcog_profiles.tsv) are written as <SMCOGxxxx>.hmm with their cutoff as GA.
PGAP models are written under their rule name (rule_suggester.pfam.libraries.pgap_rule_name), the
same name step 2 gave them, with the HMM NAME line set to it.

Writes one <domain>.hmm file per referenced profile (ASCII HMMER3 format) into <out-dir>/rule_hmms/,
matching the one-file-per-domain convention antiSMASH's subclusters module data/ directory already uses
(antismash/antismash/detection/subclusters/data/rule_hmms/) -- ready to drop straight into it --
plus <out-dir>/hmmdetails.txt for that data/ directory and <out-dir>/cutoffs.tsv (the exact cutoffs,
for review). <out-dir> thus mirrors the module's data/ layout.

The .hmm files are written as the databases ship them, except for the NAME line of PGAP and smCOG
models, which is set to the rule name (antiSMASH matches hits to hmmdetails.txt by NAME). Cutoffs go
only into hmmdetails.txt and cutoffs.tsv -- the same cutoffs step 2 used:
  - a cutoff override (data/pfam/cutoff_overrides.tsv), if the model has one
  - smCOGs: the curated cutoff (data/smcogs/smcog_profiles.tsv)
  - otherwise the model's own GA (Pfam, TIGRFAM, NCBIFAM)
hmmdetails.txt rounds the cutoff down to a whole number, since antiSMASH reads it with int(); rounding
down keeps every hit step 2 counted. Don't use antiSMASH's data/extract_details.py on these files: it
takes the TC line, which ignores overrides and is missing for smCOGs.

"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.export.parse import extract_domain_ids
from rule_suggester.pfam.cutoffs import load_cutoff_overrides
from rule_suggester.pfam.libraries import PgapLibrary, SmcogLibrary, iter_hmm_file, load_pfam_names, load_smcog_profiles


def _text(value: bytes | str) -> str:
    return value.decode() if isinstance(value, bytes) else value

DATA = Path(__file__).resolve().parent.parent / "data"
DEFAULT_PFAM_HMM = DATA / "pfam" / "Pfam-A.hmm"
DEFAULT_PFAM_CLANS = DATA / "pfam" / "Pfam-A.clans.tsv"
DEFAULT_PGAP_HMM = DATA / "ncbi_pgap" / "hmm_PGAP.HMM.tgz"
DEFAULT_PGAP_METADATA = DATA / "ncbi_pgap" / "hmm_PGAP.tsv"
DEFAULT_SMCOG_HMM = DATA / "smcogs" / "smcogs.hmm"
DEFAULT_SMCOG_PROFILES = DATA / "smcogs" / "smcog_profiles.tsv"
DEFAULT_CUTOFF_OVERRIDES = DATA / "pfam" / "cutoff_overrides.tsv"
_DESC_PREFIX = re.compile(r"^(NCBIFAM|TIGRFAM|JCVI|PGAP)\s*:\s*")


def hmmdetails_line(name: str, description: str, cutoff: float, accession: str) -> str:
    """One antiSMASH subclusters hmmdetails.txt line: name, description, whole-number cutoff (rounded
    down), path relative to data/, accession (with version; empty if the model has none)."""
    description = " ".join(_DESC_PREFIX.sub("", description).split()) or name
    return "\t".join([name, description, str(math.floor(cutoff + 1e-9)), f"rule_hmms/{name}.hmm", accession])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rules-file", type=Path, required=True,
                         help="A combined rules file, e.g. output/subcluster_rules.txt.")
    parser.add_argument("--pfam-hmm", type=Path, default=DEFAULT_PFAM_HMM)
    parser.add_argument("--pfam-clans", type=Path, default=DEFAULT_PFAM_CLANS,
                         help="Pfam-A.clans.tsv, for the Pfam names PGAP rule names must not collide with.")
    parser.add_argument("--pgap-hmm", type=Path, default=DEFAULT_PGAP_HMM)
    parser.add_argument("--pgap-metadata", type=Path, default=DEFAULT_PGAP_METADATA)
    parser.add_argument("--no-pgap", action="store_true", help="Don't look in PGAP.")
    parser.add_argument("--smcog-hmm", type=Path, default=DEFAULT_SMCOG_HMM)
    parser.add_argument("--smcog-profiles", type=Path, default=DEFAULT_SMCOG_PROFILES,
                         help="Curated smCOG list (as for step 2); its profiles are written as <SMCOGxxxx>.hmm "
                              "with the curated cutoff as GA. Used if the file exists.")
    parser.add_argument("--cutoff-overrides", type=Path, default=DEFAULT_CUTOFF_OVERRIDES,
                         help="Cutoff overrides as for step 2; used if the file exists.")
    parser.add_argument("--out-dir", type=Path, required=True,
                         help="Output directory, e.g. output/. Gets hmmdetails.txt and cutoffs.tsv, and "
                              "one <domain>.hmm file per referenced profile in its rule_hmms/ subdirectory.")
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

    overrides = load_cutoff_overrides(args.cutoff_overrides) if args.cutoff_overrides.exists() else {}
    hmm_dir = args.out_dir / "rule_hmms"
    hmm_dir.mkdir(parents=True, exist_ok=True)
    found: set[str] = set()
    cutoff_rows = []
    details_lines = []
    libraries = [("Pfam", args.pfam_hmm.name, lambda: iter_hmm_file(args.pfam_hmm))]
    if not args.no_pgap:
        pgap = PgapLibrary.load(args.pgap_hmm, args.pgap_metadata, load_pfam_names(args.pfam_clans))
        libraries.append(("PGAP", args.pgap_hmm.name, pgap.iter_hmms))
    if args.smcog_profiles.exists() and any(name.startswith("SMCOG") for name in wanted):
        smcog = SmcogLibrary(args.smcog_hmm, load_smcog_profiles(args.smcog_profiles))
        libraries.append(("smCOG", args.smcog_hmm.name, lambda: smcog.iter_hmms(apply_cutoffs=False)))
    for kind, label, hmms in libraries:
        if found == wanted:
            break
        print(f"Scanning {label} ...")
        for hmm in hmms():
            name = _text(hmm.name)
            if name in wanted and name not in found:
                acc = _text(hmm.accession or "").split(".")[0]
                ga = hmm.cutoffs.gathering if hmm.cutoffs.gathering_available() else (float("nan"),) * 2
                if acc in overrides:
                    cutoff, reason = (overrides[acc],) * 2, "cutoff override"
                elif kind == "PGAP":
                    cutoff, reason = ga, f"{pgap.model(_text(hmm.accession)).source} GA"
                elif kind == "smCOG":
                    cutoff, reason = (smcog.profiles[name],) * 2, "curated smCOG cutoff"
                else:
                    cutoff, reason = ga, "Pfam GA"
                out_path = hmm_dir / f"{name}.hmm"
                with out_path.open("wb") as fh:
                    hmm.write(fh, binary=False)
                cutoff_rows.append([name, acc, kind, f"{ga[1]:g}", f"{cutoff[1]:.2f}", reason])
                details_lines.append(hmmdetails_line(
                    name, smcog.descriptions[name] if kind == "smCOG" else _text(hmm.description or ""), cutoff[1],
                    "" if kind == "smCOG" else _text(hmm.accession or ""),
                ))
                found.add(name)
                if found == wanted:
                    break
    with (args.out_dir / "cutoffs.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["name", "accession", "library", "model_ga", "cutoff", "reason"])
        w.writerows(sorted(cutoff_rows))
    (args.out_dir / "hmmdetails.txt").write_text("\n".join(sorted(details_lines, key=str.lower)) + "\n")

    missing = sorted(wanted - found)
    print(f"{len(found)}/{len(wanted)} written to {hmm_dir}; hmmdetails.txt and cutoffs.tsv in {args.out_dir}")
    if missing:
        print(f"WARNING: {len(missing)} domain(s) not found in {' / '.join(l for _, l, _ in libraries)}:", file=sys.stderr)
        for name in missing:
            print(f"  {name}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
