#!/usr/bin/env python3
"""CLI entry point for building a full antiSMASH RULE block (Phase 3, building on
scripts/pfam_domains.py's Phase 2 output).

RULE = the SCG id, e.g. SCG0031.
CATEGORY = always "Subcluster".
DESCRIPTION = the curated product name from data/subcluster_groups.json.
CUTOFF = estimated from real gene spacing in data/genbanks/<group>/<bgc_id>.gbk, when available
(see rule_suggester.rules.gene_ranges); NEIGHBOURHOOD is always 0 (see render_rule_block()).
CONDITIONS = core groups' representative domains, cds()-ANDed per group, ANDed across groups.
EXTENDERS = non-core groups' representative domains, cds()-ANDed per group, ORed across groups.

Also writes rule_conditions.md: the same block plus the per-group detail and warnings behind it,
for review -- the domain selection (CONDITIONS/EXTENDERS content) is not guaranteed correct and
may need manual editing; that detail is there to make checking it easier.

--segments-dir accepts the same single-group-or-parent-directory semantics as the other scripts.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.phase1_io import load_groups_tsv
from rule_suggester.rules.conditions import GroupSummary, RuleConditions, build_conditions, build_group_summaries
from rule_suggester.rules.gene_ranges import BgcSpan, compute_bgc_spans, estimate_cutoff_kb
from rule_suggester.rules.product_names import load_product_names, sanitize_description

DEFAULT_GENBANK_ROOT = Path(__file__).resolve().parent.parent / "data" / "genbanks"
DEFAULT_SUBCLUSTER_GROUPS_JSON = Path(__file__).resolve().parent.parent / "data" / "subcluster_groups.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segments-dir", "--segment-dir", dest="segments_dir", type=Path, required=True,
                         help="A single group's phase-1/2 output dir (e.g. output/segments/SCG0001_streptidine), "
                              "or a parent dir of them (e.g. output/segments).")
    parser.add_argument(
        "--genbank-dir", type=Path, default=None,
        help="Single group only: dir of <bgc_id>.gbk files for this product group's example BGCs, "
             "used to estimate CUTOFF from real gene spacing. Overrides --genbank-root.",
    )
    parser.add_argument(
        "--genbank-root", type=Path, default=DEFAULT_GENBANK_ROOT,
        help="Parent of per-group GenBank dirs; each group uses <genbank-root>/<group dir name>.",
    )
    parser.add_argument(
        "--subcluster-groups-json", type=Path, default=DEFAULT_SUBCLUSTER_GROUPS_JSON,
        help="Curated product-name lookup for RULE DESCRIPTION, keyed by SCG id.",
    )
    return parser.parse_args()


def _format_conditions_multiline(terms: list[str], keyword: str) -> list[str]:
    """One top-level term per line, continuation lines aligned under where the first term
    starts: matches manually-written rules' style, and keeps long AND-chains scannable.
    """
    indent = " " * (4 + len(keyword) + 1)
    lines = [f"    {keyword} {terms[0]}"]
    lines.extend(f"{indent}and {term}" for term in terms[1:])
    return lines


def render_rule_block(rule_id: str, description: str, cutoff: int | None, rule: RuleConditions) -> list[str]:
    lines = [
        f"RULE {rule_id}",
        "    CATEGORY Subcluster",
        f"    DESCRIPTION {sanitize_description(description)}",
        f"    CUTOFF {cutoff if cutoff is not None else '<unresolved -- see CUTOFF/NEIGHBOURHOOD estimate below>'}",
        "    NEIGHBOURHOOD 0",
    ]
    if rule.condition_terms:
        lines.extend(_format_conditions_multiline(rule.condition_terms, "CONDITIONS"))
    else:
        lines.append("    CONDITIONS <none -- no core group has a representative domain>")
    if rule.extenders:
        lines.append(f"    EXTENDERS {rule.extenders}")
    return lines


def _group_row(s: GroupSummary, include_bgcs: bool) -> str:
    genes = ", ".join(s.gene_names) or "--"
    domains = " + ".join(s.representative_domains) or "--"
    contributes = "yes" if s.contributes else "no (no signature domain)"
    if include_bgcs:
        return f"| {s.group_id} | {genes} | {s.n_bgcs}/{s.n_bgcs_total} | {domains} | {contributes} |"
    return f"| {s.group_id} | {genes} | {domains} | {contributes} |"


def build_warnings(
    core: list[GroupSummary],
    extra: list[GroupSummary],
    rule: RuleConditions,
    unresolved_bgcs: list[str],
    partial_bgcs: list[BgcSpan],
    description: str | None,
    rule_id: str,
    subcluster_groups_json: Path,
) -> list[str]:
    warnings: list[str] = []

    if description is None:
        try:
            json_display = subcluster_groups_json.relative_to(Path.cwd())
        except ValueError:
            json_display = subcluster_groups_json
        warnings.append(
            f"**No product name found for {rule_id} in `{json_display}`** "
            "-- DESCRIPTION left as a placeholder."
        )
    if unresolved_bgcs:
        warnings.append(f"**No GenBank file found for BGC(s):** {', '.join(unresolved_bgcs)}")
    for s in partial_bgcs:
        warnings.append(f"**Genes not found in {s.bgc_id}'s GenBank file:** {', '.join(s.unmatched)}")
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
    return warnings


def render_rule_block_with_warnings(
    rule_id: str, description: str, cutoff: int | None, rule: RuleConditions, warnings: list[str]
) -> list[str]:
    """The rule.txt content: the RULE block, preceded by a '#'-commented summary of any
    warnings, so a rule copied out of this file on its own still carries that context
    instead of it being stranded in rule_conditions.md.
    """
    lines: list[str] = []
    if warnings:
        lines.append(f"# WARNING: {len(warnings)} warning(s) below -- see rule_conditions.md for full detail.")
        lines.extend(f"# - {w.replace('**', '')}" for w in warnings)
        lines.append("")
    lines.extend(render_rule_block(rule_id, description, cutoff, rule))
    return lines


def render(
    group_name: str,
    summaries: list[GroupSummary],
    rule: RuleConditions,
    spans: list[BgcSpan],
    cutoff: int | None,
    genbank_dir: Path,
    rule_id: str,
    description: str | None,
    warnings: list[str],
) -> str:
    core = [s for s in summaries if s.is_core]
    extra = [s for s in summaries if not s.is_core]

    lines = [f"# Rule conditions: {group_name}", ""]

    lines.append("## RULE block")
    lines.append("")
    lines.append(
        "Ready to paste into an antiSMASH rule file. CONDITIONS/EXTENDERS content is built "
        "automatically and is not guaranteed correct -- review against the detail and warnings "
        "below before using it."
    )
    lines.append("")
    lines.append("```")
    lines.extend(render_rule_block(rule_id, description or "<unknown -- not in subcluster_groups.json>", cutoff, rule))
    lines.append("```")
    lines.append("")

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

    lines.append("## CUTOFF / NEIGHBOURHOOD estimate")
    lines.append("")
    try:
        genbank_dir_display = genbank_dir.relative_to(Path.cwd())
    except ValueError:
        genbank_dir_display = genbank_dir
    if spans:
        lines.append(
            "antiSMASH looks CUTOFF kb to the left and right of each anchor CDS, so CUTOFF is a radius, "
            "not the width of the whole locus. Per reference BGC in "
            f"`{genbank_dir_display}`, the radius is the smallest reach (edge to edge) that one "
            "CONDITIONS gene needs to cover every other CONDITIONS gene. CUTOFF is the largest radius "
            "across BGCs, rounded up to the nearest 5 kb -- below that, "
            "some example's CONDITIONS genes never all fall within CUTOFF of a single anchor. The "
            "full span (all genes in the group, core and non-core) is shown for reference only. "
            "NEIGHBOURHOOD is always 0: CONDITIONS/EXTENDERS already require every group's own "
            "domains, so there's no extra accessory search radius left to justify."
        )
        lines.append("")
        lines.append("| BGC | genes | matched in GenBank | span (kb) | radius (kb) |")
        lines.append("|---|---|---|---|---|")
        for s in spans:
            span_str = f"{s.span_kb:.1f}" if s.span_kb is not None else "unresolved"
            radius_str = f"{s.radius_kb:.1f}" if s.radius_kb is not None else "unresolved"
            lines.append(f"| {s.bgc_id} | {s.n_genes} | {s.n_matched}/{s.n_genes} | {span_str} | {radius_str} |")
        lines.append("")
        if cutoff is not None:
            lines.append("```")
            lines.append(f"CUTOFF {cutoff}")
            lines.append("NEIGHBOURHOOD 0")
            lines.append("```")
        else:
            lines.append(
                "**WARNING: no BGC's genes could be resolved against its GenBank file -- "
                "CUTOFF can't be estimated.**"
            )
    else:
        lines.append(
            f"*No GenBank directory found at `{genbank_dir_display}` -- CUTOFF can't be estimated. "
            f"Pass --genbank-dir, or add data/genbanks/{group_name}/<bgc_id>.gbk files.*"
        )
    lines.append("")

    lines.append("## Warnings")
    lines.append("")
    lines.extend(warnings if warnings else ["None."])
    lines.append("")

    return "\n".join(lines)


def _is_group_dir(path: Path) -> bool:
    return (path / "groups.tsv").exists()


def discover_group_dirs(segments_dir: Path) -> list[Path]:
    if not segments_dir.is_dir():
        raise ValueError(f"--segments-dir {segments_dir} is not a directory")
    if _is_group_dir(segments_dir):
        return [segments_dir]
    group_dirs = sorted(p for p in segments_dir.iterdir() if p.is_dir() and _is_group_dir(p))
    if not group_dirs:
        raise ValueError(f"{segments_dir} contains no group directories with a groups.tsv")
    return group_dirs


def build_group_rule(
    segment_dir: Path, genbank_dir: Path, subcluster_groups_json: Path, product_names: dict[str, str]
) -> tuple[str, list[str]]:
    """Build one group's rule_conditions.md and rule.txt; returns (rule_conditions text, warnings).
    Raises FileNotFoundError if the group's phase-1/2 inputs are incomplete.
    """
    groups_tsv = segment_dir / "groups.tsv"
    presence_tsv = segment_dir / "presence.tsv"
    group_domains_tsv = segment_dir / "group_domains.tsv"
    missing = [p.name for p in (groups_tsv, presence_tsv, group_domains_tsv) if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"{segment_dir} is missing {', '.join(missing)} (run segment_groups.py and pfam_domains.py first)"
        )

    summaries = build_group_summaries(groups_tsv, presence_tsv, group_domains_tsv)
    rule = build_conditions(summaries)

    spans: list[BgcSpan] = []
    cutoff = None
    if genbank_dir.is_dir():
        memberships = load_groups_tsv(groups_tsv)
        condition_groups = {s.group_id for s in summaries if s.is_core and s.contributes}
        spans = compute_bgc_spans(memberships, genbank_dir, condition_groups)
        cutoff = estimate_cutoff_kb(spans)

    rule_id = segment_dir.name.split("_", 1)[0]
    description = product_names.get(rule_id)

    core = [s for s in summaries if s.is_core]
    extra = [s for s in summaries if not s.is_core]
    unresolved_bgcs = [s.bgc_id for s in spans if s.span_kb is None]
    partial_bgcs = [s for s in spans if s.span_kb is not None and s.unmatched]
    warnings = build_warnings(
        core, extra, rule, unresolved_bgcs, partial_bgcs, description, rule_id, subcluster_groups_json,
    )

    text = render(
        segment_dir.name, summaries, rule, spans, cutoff, genbank_dir,
        rule_id, description, warnings,
    )
    (segment_dir / "rule_conditions.md").write_text(text)

    rule_block_lines = render_rule_block_with_warnings(
        rule_id, description or "<unknown>", cutoff, rule, warnings,
    )
    (segment_dir / "rule.txt").write_text("\n".join(rule_block_lines) + "\n")
    return text, warnings


def main() -> int:
    args = parse_args()
    try:
        group_dirs = discover_group_dirs(args.segments_dir)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    single_group = len(group_dirs) == 1 and group_dirs[0] == args.segments_dir
    if args.genbank_dir is not None and not single_group:
        print("error: --genbank-dir only applies to a single group; use --genbank-root instead", file=sys.stderr)
        return 1

    product_names = load_product_names(args.subcluster_groups_json) if args.subcluster_groups_json.exists() else {}

    failed = []
    for group_dir in group_dirs:
        genbank_dir = args.genbank_dir or (args.genbank_root / group_dir.name)
        try:
            text, warnings = build_group_rule(group_dir, genbank_dir, args.subcluster_groups_json, product_names)
        except FileNotFoundError as e:
            print(f"error: {e}", file=sys.stderr)
            failed.append(group_dir.name)
            continue

        if single_group:
            print(text)
            print(f"Written to {group_dir / 'rule_conditions.md'} and {group_dir / 'rule.txt'}", file=sys.stderr)
        else:
            status = f"{len(warnings)} warning(s)" if warnings else "ok"
            print(f"  {group_dir.name}: {status}")

    if not single_group:
        print(f"Done: {len(group_dirs) - len(failed)} group(s) built, {len(failed)} failed.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
