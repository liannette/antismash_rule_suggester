#!/usr/bin/env python3
"""CLI entry point for building compound_details.txt (antiSMASH subclusters module format) from
a combined rules file (scripts/combine_rules.py's output) and the curated product metadata in
data/subcluster_groups.json.

Format (see antismash/antismash/detection/subclusters/compounds.py): tab-separated, one row per
RULE id with columns rule_name, compound_name, smiles, classification (semicolon-separated).
Every RULE in the target rules file needs a matching row, or antiSMASH's subclusters module
KeyErrors loading it -- this writes exactly one row per RULE id found in --rules-file, in the
order they appear there.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.export.compounds import load_products
from rule_suggester.export.parse import extract_rule_ids

DEFAULT_SUBCLUSTER_GROUPS_JSON = Path(__file__).resolve().parent.parent / "data" / "subcluster_groups.json"
DEFAULT_OUT = Path(__file__).resolve().parent.parent / "output" / "compound_details.txt"

HEADER = "# rule_name\tcompound_name\tsmiles\tclassification (semicolon-separated)"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rules-file", type=Path, required=True,
                         help="A combined rules file, e.g. output/subcluster_rules.txt.")
    parser.add_argument("--subcluster-groups-json", type=Path, default=DEFAULT_SUBCLUSTER_GROUPS_JSON)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.rules_file.exists():
        print(f"error: {args.rules_file} not found", file=sys.stderr)
        return 1
    if not args.subcluster_groups_json.exists():
        print(f"error: {args.subcluster_groups_json} not found", file=sys.stderr)
        return 1

    rule_ids = extract_rule_ids(args.rules_file.read_text())
    products = load_products(args.subcluster_groups_json)

    lines = [HEADER]
    missing: list[str] = []
    incomplete: list[str] = []
    for rule_id in rule_ids:
        product = products.get(rule_id)
        if product is None:
            missing.append(rule_id)
            continue
        name = product.get("name") or ""
        smiles = product.get("smiles") or ""
        classification = ";".join(product.get("classification") or [])
        if not smiles or not classification:
            incomplete.append(rule_id)
        lines.append(f"{rule_id}\t{name}\t{smiles}\t{classification}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")

    written = len(rule_ids) - len(missing)
    print(f"{written}/{len(rule_ids)} rule(s) written to {args.out}")
    if incomplete:
        print(f"WARNING: {len(incomplete)} rule(s) have an empty smiles or classification:", file=sys.stderr)
        for rule_id in incomplete:
            print(f"  {rule_id}", file=sys.stderr)
    if missing:
        print(
            f"WARNING: {len(missing)} rule(s) have no entry in {args.subcluster_groups_json} "
            "-- compound_details.txt is incomplete, antiSMASH will KeyError loading these rules:",
            file=sys.stderr,
        )
        for rule_id in missing:
            print(f"  {rule_id}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
