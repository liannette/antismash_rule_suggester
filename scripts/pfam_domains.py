#!/usr/bin/env python3
"""CLI entry point for Pfam/PGAP domain mapping (docs/segment_grouping_spec.md, Phase 2).

Reads phase-1 outputs (segments.tsv / groups.tsv, from scripts/segment_groups.py) and maps domains
-- found by hmmsearch of Pfam and the NCBI PGAP library (TIGRFAM, NCBIFAM, PRK) on the original
full-length proteins -- onto them. PGAP is searched unless --no-pgap is given. --fasta-dir has
the same single-group-or-parent-directory semantics as scripts/segment_groups.py.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.pfam.cutoffs import load_cutoff_overrides
from rule_suggester.pfam.libraries import PgapLibrary
from rule_suggester.pfam.params import PfamParams
from rule_suggester.pfam.pipeline import run_pfam

DATA_PFAM_DIR = Path(__file__).resolve().parent.parent / "data" / "pfam"
DEFAULT_PFAM_HMM = DATA_PFAM_DIR / "Pfam-A.hmm"
DEFAULT_PFAM_CLANS = DATA_PFAM_DIR / "Pfam-A.clans.tsv"
DEFAULT_CUTOFF_OVERRIDES = DATA_PFAM_DIR / "cutoff_overrides.tsv"
DATA_PGAP_DIR = Path(__file__).resolve().parent.parent / "data" / "ncbi_pgap"
DEFAULT_PGAP_HMM = DATA_PGAP_DIR / "hmm_PGAP.HMM.tgz"
DEFAULT_PGAP_METADATA = DATA_PGAP_DIR / "hmm_PGAP.tsv"


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
    parser.add_argument("--pgap-hmm", type=Path, default=DEFAULT_PGAP_HMM,
                         help="NCBI PGAP HMM library: the hmm_PGAP.HMM.tgz archive (read without unpacking) "
                              "or a concatenated HMMER3 file. Default: %(default)s")
    parser.add_argument("--pgap-metadata", type=Path, default=DEFAULT_PGAP_METADATA,
                         help="hmm_PGAP.tsv (labels, family types, product names). Default: %(default)s")
    parser.add_argument("--no-pgap", action="store_true", help="Search Pfam only.")
    parser.add_argument("--cutoff-overrides", type=Path, default=DEFAULT_CUTOFF_OVERRIDES,
                         help="TSV (pfam_acc, cutoff, note) of Pfam/PGAP models whose GA is replaced by a bitscore "
                              "cutoff. Used if the file exists; default: %(default)s")
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

    if not args.no_pgap:
        for path in (args.pgap_hmm, args.pgap_metadata):
            if not path.exists():
                print(f"PGAP file not found at {path}. Stopping (use --no-pgap to search Pfam only).",
                      file=sys.stderr)
                return 1

    overrides: dict[str, float] = {}
    if args.cutoff_overrides.exists():
        try:
            overrides = load_cutoff_overrides(args.cutoff_overrides)
        except ValueError as e:
            print(f"error: {e}", file=sys.stderr)
            return 1
        known = {line.split("\t", 1)[0] for line in args.pfam_clans.open() if line.strip()}
        if not args.no_pgap:
            known |= PgapLibrary.load(args.pgap_hmm, args.pgap_metadata, set()).accessions()
        unknown = sorted(set(overrides) - known)
        if unknown:
            print(f"error: {args.cutoff_overrides} lists model(s) in neither Pfam nor PGAP: "
                  f"{', '.join(unknown)}", file=sys.stderr)
            return 1
        print(f"Cutoff overrides from {args.cutoff_overrides}: "
              + (", ".join(f"{a} {c:g}" for a, c in sorted(overrides.items())) or "none"))
    else:
        print(f"No cutoff overrides file at {args.cutoff_overrides}: every Pfam uses GA.")

    params = PfamParams(
        hmm_path=args.pfam_hmm,
        min_prevalence=args.min_prevalence,
        straddle_threshold=args.straddle_threshold,
        cpus=args.cpus,
        cutoff_overrides=overrides,
        pgap_hmm_path=None if args.no_pgap else args.pgap_hmm,
        pgap_metadata_path=None if args.no_pgap else args.pgap_metadata,
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
