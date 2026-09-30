#!/usr/bin/env python3
"""CLI entry point for building antiSMASH-style CONDITIONS/EXTENDERS fragments from one segment
group's representative Pfam domains (Phase 3, building on scripts/pfam_domains.py's Phase 2
output).

CONDITIONS = core groups' representative domains, cds()-ANDed per group, ANDed across groups.
EXTENDERS = non-core groups' representative domains, cds()-ANDed per group, ORed across groups.

Writes only these two fragments plus the per-group detail behind them -- not a complete RULE
block, since CATEGORY/CUTOFF/NEIGHBOURHOOD aren't derivable from segment/domain data and a
placeholder default would misrepresent them as decided; fill those in by hand when wiring the
fragments into an actual antiSMASH rule file.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.rules.conditions import GroupSummary, RuleConditions, build_conditions, build_group_summaries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segment-dir", type=Path, required=True,
                         help="A single group's phase-1/2 output dir, e.g. output/segments/SCG0001_streptidine.")
    return parser.parse_args()


def _group_row(s: GroupSummary, include_bgcs: bool) -> str:
    genes = ", ".join(s.gene_names) or "--"
    domains = " + ".join(s.representative_domains) or "--"
    contributes = "yes" if s.contributes else "no (no signature domain)"
    if include_bgcs:
        return f"| {s.group_id} | {genes} | {s.n_bgcs}/{s.n_bgcs_total} | {domains} | {contributes} |"
    return f"| {s.group_id} | {genes} | {domains} | {contributes} |"


def render(group_name: str, summaries: list[GroupSummary], rule: RuleConditions) -> str:
    core = [s for s in summaries if s.is_core]
    extra = [s for s in summaries if not s.is_core]

    lines = [f"# Rule conditions: {group_name}", ""]

    lines.append("## Core groups (CONDITIONS)")
    lines.append("")
    lines.append(
        "Every core group must contribute a hit. A group with no signature domain in Step 2 "
        "(no domain at 100% prevalence) has nothing to contribute and is left out of CONDITIONS "
        "below -- that silently weakens the rule, so check the 'no' rows (also listed under "
        "Warnings)."
    )
    lines.append("")
    if core:
        lines.append("| Group | Genes | Representative domains | Contributes |")
        lines.append("|---|---|---|---|")
        lines.extend(_group_row(s, include_bgcs=False) for s in core)
    else:
        lines.append("None.")
    lines.append("")
    if rule.repeated_domains:
        for domain, group_ids in rule.repeated_domains.items():
            lines.append(
                f"Note: {domain} is the sole representative in {len(group_ids)} core groups "
                f"({', '.join(group_ids)}) -- folded into `repeatable({len(group_ids)}, [{domain}])` "
                "to require that many separate genes, instead of one deduped 'and' term."
            )
        lines.append("")
    if rule.conditions:
        lines.append("```")
        lines.append(f"CONDITIONS {rule.conditions}")
        lines.append("```")
    else:
        lines.append("**WARNING: no core group has a representative domain -- CONDITIONS can't be built.**")
    lines.append("")

    lines.append("## Non-core groups (EXTENDERS)")
    lines.append("")
    if extra:
        lines.append("| Group | Genes | BGCs | Representative domains | Contributes |")
        lines.append("|---|---|---|---|---|")
        lines.extend(_group_row(s, include_bgcs=True) for s in extra)
        lines.append("")
        if rule.extenders:
            lines.append("```")
            lines.append(f"EXTENDERS {rule.extenders}")
            lines.append("```")
        else:
            lines.append("**WARNING: no non-core group has a representative domain -- EXTENDERS can't be built.**")
    else:
        lines.append("None -- every group is core.")
    lines.append("")

    lines.append("## Warnings")
    lines.append("")
    warnings: list[str] = []

    missing_core = [s.group_id for s in core if not s.contributes]
    if missing_core:
        warnings.append(
            f"**No signature domain, excluded from CONDITIONS:** {', '.join(missing_core)}"
        )
    missing_extra = [s.group_id for s in extra if not s.contributes]
    if missing_extra:
        warnings.append(
            f"**No signature domain, excluded from EXTENDERS:** {', '.join(missing_extra)}"
        )
    if not core:
        warnings.append("**No core groups at all** -- CONDITIONS can't anchor a rule on this product group.")
    if rule.duplicate_condition_groups:
        for term, group_ids in rule.duplicate_condition_groups.items():
            warnings.append(
                f"**Duplicate CONDITIONS term across different core groups:** {', '.join(group_ids)} "
                f"all reduce to `{term}` -- can't be folded into repeatable() since it has more than "
                "one domain; check whether these are genuinely distinct genes."
            )
    if rule.redundant_extenders:
        for dominated, dominating in rule.redundant_extenders:
            warnings.append(
                f"**Redundant EXTENDERS clause:** {dominated}'s representative domains are a superset "
                f"of {dominating}'s, so {dominated}'s term can never independently satisfy EXTENDERS "
                f"-- {dominating} alone already covers it."
            )
    lines.extend(warnings if warnings else ["None."])
    lines.append("")

    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    groups_tsv = args.segment_dir / "groups.tsv"
    presence_tsv = args.segment_dir / "presence.tsv"
    group_domains_tsv = args.segment_dir / "group_domains.tsv"
    missing = [p.name for p in (groups_tsv, presence_tsv, group_domains_tsv) if not p.exists()]
    if missing:
        print(
            f"error: {args.segment_dir} is missing {', '.join(missing)} "
            "(run segment_groups.py and pfam_domains.py first)",
            file=sys.stderr,
        )
        return 1

    summaries = build_group_summaries(groups_tsv, presence_tsv, group_domains_tsv)
    rule = build_conditions(summaries)
    text = render(args.segment_dir.name, summaries, rule)

    out_path = args.segment_dir / "rule_conditions.md"
    out_path.write_text(text)
    print(text)
    print(f"Written to {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
