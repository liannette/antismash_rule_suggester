#!/usr/bin/env python3
"""CLI entry point for building a full antiSMASH RULE block (Phase 3, building on
scripts/pfam_domains.py's Phase 2 output).

RULE = the bb_id (building block), e.g. BB0003.
CATEGORY = always "Subcluster".
DESCRIPTION = the curated building-block name from data/building_blocks.json.
CUTOFF = estimated from real gene spacing in data/genbanks/<group>/<sc_id>.gbk, when available
(see rule_suggester.rules.gene_ranges); NEIGHBOURHOOD is always 0 (see render_rule_block()).
CONDITIONS = core groups' representative domains, cds()-ANDed per group, ANDed across groups.
EXTENDERS = non-core groups' representative domains, cds()-ANDed per group, ORed across groups.

Also writes rule_conditions.md: the same block plus the per-group detail and warnings behind it,
for review -- the domain selection (CONDITIONS/EXTENDERS content) is not guaranteed correct and
may need manual editing; that detail is there to make checking it easier.

--segments-dir accepts the same single-group-or-parent-directory semantics as the other scripts.
"""

import argparse
import csv
import functools
import re
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.phase1_io import load_groups_tsv
from rule_suggester.rules.conditions import GroupSummary, RuleConditions, build_conditions, build_group_summaries
from rule_suggester.export.parse import _IDENTIFIER_RE, _KEYWORDS
from rule_suggester.rules.gene_ranges import ScSpan, compute_sc_spans, estimate_cutoff_kb
from rule_suggester.rules.manual import ManualRule, load_manual_rules
from rule_suggester.rules.product_names import load_product_names, sanitize_description

DEFAULT_GENBANK_ROOT = Path(__file__).resolve().parent.parent / "data" / "genbanks"
DEFAULT_BUILDING_BLOCKS_JSON = Path(__file__).resolve().parent.parent / "data" / "building_blocks.json"
DEFAULT_MANUAL_RULES = Path(__file__).resolve().parent.parent / "data" / "manual_rules.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segments-dir", "--segment-dir", dest="segments_dir", type=Path, required=True,
                         help="A single group's phase-1/2 output dir (e.g. output/segments/BB0003_L-dihydrostreptose), "
                              "or a parent dir of them (e.g. output/segments).")
    parser.add_argument(
        "--genbank-dir", type=Path, default=None,
        help="Single group only: dir of <sc_id>.gbk files for this product group's example subclusters, "
             "used to estimate CUTOFF from real gene spacing. Overrides --genbank-root.",
    )
    parser.add_argument(
        "--genbank-root", type=Path, default=DEFAULT_GENBANK_ROOT,
        help="Parent of per-group GenBank dirs; each group uses <genbank-root>/<group dir name>.",
    )
    parser.add_argument(
        "--building-blocks-json", type=Path, default=DEFAULT_BUILDING_BLOCKS_JSON,
        help="Curated building-block names for RULE DESCRIPTION, keyed by bb_id (default: %(default)s).",
    )
    parser.add_argument(
        "--manual-rules", type=Path, default=DEFAULT_MANUAL_RULES,
        help="Hand-written CONDITIONS/EXTENDERS/CUTOFF/DROP per rule (rule_suggester/rules/manual.py). "
             "Used if the file exists; default: %(default)s",
    )
    return parser.parse_args()


def _identifiers(condition: str) -> set[str]:
    return set(_IDENTIFIER_RE.findall(condition)) - _KEYWORDS


_NEGATED_RE = re.compile(r"\bnot\s+([A-Za-z0-9_-]+)")


def _negated_identifiers(condition: str) -> set[str]:
    """Identifiers directly after `not` (e.g. `not SMCOG1147`, `cds(AMP-binding and not ketoacyl-synt)`):
    models that must be absent, so they needn't be a signature of this building block's own genes."""
    return set(_NEGATED_RE.findall(condition)) - _KEYWORDS


@functools.cache
def searched_models(segments_root: Path) -> frozenset[str]:
    """Every model step 2 searched with: those it reported a hit for in any building block under
    segments_root, plus the curated smCOG and NRPS/PKS profiles it was configured to search.

    A curated profile is searched whether or not it hits anything, and a cutoff chosen to exclude
    something (e.g. Condensation_LCL at 150, which no reference A-PCP reaches) is *meant* to leave
    no hits. Such a profile is still a valid `not` term, so the configured lists count as searched;
    otherwise a correct exclusion is reported as a typo.
    """
    names: set[str] = set()
    for domains_tsv in segments_root.glob("*/domains.tsv"):
        with domains_tsv.open(newline="") as f:
            names.update(row["pfam_name"] for row in csv.DictReader(f, delimiter="\t"))
    data = Path(__file__).resolve().parent.parent / "data"
    for profile_list, column in ((data / "smcogs" / "smcog_profiles.tsv", "smcog_id"),
                                 (data / "nrps_pks" / "nrps_pks_profiles.tsv", "profile")):
        if not profile_list.exists():
            continue
        with profile_list.open(newline="") as f:
            rows = csv.DictReader((line for line in f if line.strip() and not line.startswith("#")),
                                  delimiter="\t")
            names.update(row[column].strip() for row in rows if row.get(column))
    return frozenset(names)


def _format_conditions_multiline(terms: list[str], keyword: str) -> list[str]:
    """One top-level term per line, continuation lines aligned under where the first term
    starts: matches manually-written rules' style, and keeps long AND-chains scannable.
    """
    indent = " " * (4 + len(keyword) + 1)
    lines = [f"    {keyword} {terms[0]}"]
    lines.extend(f"{indent}and {term}" for term in terms[1:])
    return lines


def render_rule_block(rule_id: str, description: str, cutoff: int | None, rule: RuleConditions,
                      superiors: tuple[str, ...] = ()) -> list[str]:
    lines = [
        f"RULE {rule_id}",
        "    CATEGORY Subcluster",
        f"    DESCRIPTION {sanitize_description(description)}",
    ]
    # the grammar puts SUPERIORS after DESCRIPTION and before CUTOFF
    if superiors:
        lines.append(f"    SUPERIORS {', '.join(superiors)}")
    lines += [
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


def _group_row(s: GroupSummary, include_scs: bool) -> str:
    genes = ", ".join(s.gene_names) or "--"
    domains = " + ".join(s.representative_domains) or "--"
    contributes = "yes" if s.contributes else "no (no signature domain)"
    if include_scs:
        return f"| {s.group_id} | {genes} | {s.n_scs}/{s.n_scs_total} | {domains} | {contributes} |"
    return f"| {s.group_id} | {genes} | {domains} | {contributes} |"


def non_discriminating_warnings(summaries: list[GroupSummary]) -> list[str]:
    """Groups whose representative domains are all signatures of another group too.

    Representatives are picked per group, so nothing stops two groups from being represented by
    models the other group's gene also carries -- the rule then looks like it needs both genes
    while one can satisfy both terms (e.g. a dedicated epimerase whose Epimerase model the
    4,6-dehydratase gene matches as well). Reported with the alternatives, if any, that no other
    group shares, so the fix can be written by hand in data/manual_rules.tsv.
    """
    warnings: list[str] = []
    for s in summaries:
        reps = set(s.representative_domains)
        if not reps:
            continue
        covering = sorted(o.group_id for o in summaries
                          if o.group_id != s.group_id and reps <= set(o.all_signature_domains)
                          and reps != set(o.representative_domains))
        if not covering:
            continue
        shared = set().union(*(set(o.all_signature_domains) for o in summaries if o.group_id != s.group_id))
        unique = sorted(set(s.all_signature_domains) - shared)
        alternatives = (f"models only this group has: {', '.join(unique)}" if unique
                        else "no model of this group is unique to it")
        warnings.append(
            f"**Representative domain(s) of {s.group_id} are not specific to it:** "
            f"`{' + '.join(s.representative_domains)}` is also a signature of "
            f"{', '.join(covering)}, so one gene can satisfy both terms -- {alternatives}."
        )
    return warnings


# Warnings that describe the models available rather than a defect in the rule: a group whose
# representative isn't specific to it usually has no better model to choose (see the
# non_discriminating_warnings docstring), so the rule is as good as the profiles allow. Reported
# as advisory notes instead of blocking combine_rules.py.
ADVISORY_PREFIXES = ("**Representative domain(s) of",)


def split_advisory(warnings: list[str]) -> tuple[list[str], list[str]]:
    """(blocking, advisory) -- advisory ones still get reported, but don't hold the rule back."""
    advisory = [w for w in warnings if w.startswith(ADVISORY_PREFIXES)]
    return [w for w in warnings if w not in advisory], advisory


def build_warnings(
    core: list[GroupSummary],
    extra: list[GroupSummary],
    rule: RuleConditions,
    unresolved_scs: list[str],
    partial_scs: list[ScSpan],
    description: str | None,
    rule_id: str,
    building_blocks_json: Path,
    manual: ManualRule | None = None,
    known_domains: set[str] | None = None,
    searched: frozenset[str] | None = None,
) -> list[str]:
    warnings: list[str] = []
    manual_conditions = manual is not None and manual.conditions is not None

    if description is None:
        try:
            json_display = building_blocks_json.relative_to(Path.cwd())
        except ValueError:
            json_display = building_blocks_json
        warnings.append(
            f"**No product name found for {rule_id} in `{json_display}`** "
            "-- DESCRIPTION left as a placeholder."
        )
    if unresolved_scs:
        warnings.append(f"**No GenBank file found for subcluster(s):** {', '.join(unresolved_scs)}")
    for s in partial_scs:
        warnings.append(f"**Genes not found in {s.sc_id}'s GenBank file:** {', '.join(s.unmatched)}")
    missing_core = [s.group_id for s in core if not s.contributes]
    if missing_core and not manual_conditions:
        warnings.append(
            f"**No signature domain, excluded from CONDITIONS:** {', '.join(missing_core)}"
        )
    missing_extra = [s.group_id for s in extra if not s.contributes]
    if missing_extra:
        warnings.append(
            f"**No signature domain, excluded from EXTENDERS:** {', '.join(missing_extra)}"
        )
    if not core and not manual_conditions:
        warnings.append("**No core groups at all** -- CONDITIONS can't anchor a rule on this product group.")
    if not manual_conditions:
        warnings.extend(non_discriminating_warnings(core + extra))
    if rule.duplicate_condition_groups and not manual_conditions:
        for term, group_ids in rule.duplicate_condition_groups.items():
            warnings.append(
                f"**Duplicate CONDITIONS term across different core groups:** {', '.join(group_ids)} "
                f"all reduce to `{term}` -- one gene satisfies it, and it can't be folded into repeatable() "
                "automatically since it has more than one domain; write the CONDITIONS by hand in "
                "data/manual_rules.tsv if both genes should be required."
            )
    if manual is not None and known_domains is not None:
        for field_name, value in (("CONDITIONS", manual.conditions), ("EXTENDERS", manual.extenders)):
            if not value or value == "-":
                continue
            negated = _negated_identifiers(value)
            unknown = sorted(_identifiers(value) - negated - known_domains)
            if unknown:
                warnings.append(
                    f"**Manual {field_name} uses identifier(s) that are no signature model of any group here:** "
                    f"{', '.join(unknown)} -- check for typos (names as in group_domains.tsv)."
                )
            # an excluded model only has to exist: it was searched in step 2 and hit somewhere
            unknown_negated = sorted(negated - known_domains - (searched or frozenset()))
            if unknown_negated:
                warnings.append(
                    f"**Manual {field_name} excludes (`not`) model(s) step 2 never reported in any building block:** "
                    f"{', '.join(unknown_negated)} -- check for typos (names as in domains.tsv)."
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
    rule_id: str, description: str, cutoff: int | None, rule: RuleConditions, warnings: list[str],
    manual: ManualRule | None = None,
) -> list[str]:
    """The rule.txt content: the RULE block, preceded by a '#'-commented summary of any
    warnings, so a rule copied out of this file on its own still carries that context
    instead of it being stranded in rule_conditions.md.
    """
    lines: list[str] = []
    blocking, advisory = split_advisory(warnings)
    if blocking:
        lines.append(f"# WARNING: {len(blocking)} warning(s) below -- see rule_conditions.md for full detail.")
        lines.extend(f"# - {w.replace('**', '')}" for w in blocking)
        lines.append("")
    if advisory:
        lines.append(f"# NOTE: {len(advisory)} advisory note(s) -- see rule_conditions.md; "
                     "these don't hold the rule back.")
        lines.extend(f"# - {w.replace('**', '')}" for w in advisory)
        lines.append("")
    if manual is not None and not manual.drop:
        lines.append("# MANUAL (data/manual_rules.tsv): " + ("; ".join(manual.notes) or "hand-written fields"))
    lines.extend(render_rule_block(rule_id, description, cutoff, rule,
                                   manual.superiors if manual is not None else ()))
    if manual is not None and manual.drop:
        reason = "; ".join(manual.notes) or "no reason given"
        lines = [f"# DROPPED (data/manual_rules.tsv): {reason}"] + [f"# {line}" for line in lines]
    return lines


def render(
    group_name: str,
    summaries: list[GroupSummary],
    rule: RuleConditions,
    spans: list[ScSpan],
    cutoff: int | None,
    genbank_dir: Path,
    rule_id: str,
    description: str | None,
    warnings: list[str],
    manual: ManualRule | None = None,
) -> str:
    core = [s for s in summaries if s.is_core]
    extra = [s for s in summaries if not s.is_core]

    lines = [f"# Rule conditions: {group_name}", ""]
    if manual is not None:
        what = "DROPPED -- left out of combine_rules.py" if manual.drop else "hand-written fields applied"
        lines += [f"**Manual rule (data/manual_rules.tsv): {what}.**", ""]
        lines += [f"- {n}" for n in manual.notes] + [""]

    lines.append("## RULE block")
    lines.append("")
    lines.append(
        "Ready to paste into an antiSMASH rule file. CONDITIONS/EXTENDERS content is built "
        "automatically and is not guaranteed correct -- review against the detail and warnings "
        "below before using it."
    )
    lines.append("")
    lines.append("```")
    lines.extend(render_rule_block(rule_id, description or "<unknown -- not in building_blocks.json>", cutoff, rule,
                                   manual.superiors if manual is not None else ()))
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
        lines.extend(_group_row(s, include_scs=False) for s in core)
    else:
        lines.append("None.")
    lines.append("")
    merged = [s for s in core if s.member_groups]
    for s in merged:
        lines.append(
            f"Note: {', '.join(s.member_groups)} never share a subcluster but together cover all of them and share "
            f"{', '.join(sorted(s.signature_domains))} -- merged into core group {s.group_id}."
        )
    if rule.repeated_domains:
        for domain, group_ids in rule.repeated_domains.items():
            lines.append(
                f"Note: {domain} is the sole representative in {len(group_ids)} core groups "
                f"({', '.join(group_ids)}) -- folded into `repeatable({len(group_ids)}, [{domain}])` "
                "to require that many separate genes, instead of one deduped 'and' term."
            )
    if merged or rule.repeated_domains:
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
        lines.append("| Group | Genes | subclusters | Representative domains | Contributes |")
        lines.append("|---|---|---|---|---|")
        lines.extend(_group_row(s, include_scs=True) for s in extra)
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
            "not the width of the whole locus. Per reference subcluster in "
            f"`{genbank_dir_display}`, the radius is the smallest reach (edge to edge) that one "
            "CONDITIONS gene needs to cover every other CONDITIONS gene. CUTOFF is the largest radius "
            "across subclusters, rounded up to the nearest 5 kb -- below that, "
            "some example's CONDITIONS genes never all fall within CUTOFF of a single anchor. The "
            "full span (all genes in the group, core and non-core) is shown for reference only. "
            "NEIGHBOURHOOD is always 0: CONDITIONS/EXTENDERS already require every group's own "
            "domains, so there's no extra accessory search radius left to justify."
        )
        lines.append("")
        lines.append("| subcluster | genes | matched in GenBank | span (kb) | radius (kb) |")
        lines.append("|---|---|---|---|---|")
        for s in spans:
            span_str = f"{s.span_kb:.1f}" if s.span_kb is not None else "unresolved"
            radius_str = f"{s.radius_kb:.1f}" if s.radius_kb is not None else "unresolved"
            lines.append(f"| {s.sc_id} | {s.n_genes} | {s.n_matched}/{s.n_genes} | {span_str} | {radius_str} |")
        lines.append("")
        if cutoff is not None:
            lines.append("```")
            lines.append(f"CUTOFF {cutoff}")
            lines.append("NEIGHBOURHOOD 0")
            lines.append("```")
        else:
            lines.append(
                "**WARNING: no subcluster's genes could be resolved against its GenBank file -- "
                "CUTOFF can't be estimated.**"
            )
    else:
        lines.append(
            f"*No GenBank directory found at `{genbank_dir_display}` -- CUTOFF can't be estimated. "
            f"Pass --genbank-dir, or add data/genbanks/{group_name}/<sc_id>.gbk files.*"
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
    segment_dir: Path, genbank_dir: Path, building_blocks_json: Path, product_names: dict[str, str],
    manual_rules: dict[str, ManualRule] | None = None,
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

    rule_id = segment_dir.name.split("_", 1)[0]
    manual = (manual_rules or {}).get(rule_id)
    # hand-written CONDITIONS replace whatever a merge would contribute, so keep the groups as they are
    summaries = build_group_summaries(groups_tsv, presence_tsv, group_domains_tsv,
                                      merge_complementary=manual is None or manual.conditions is None)
    rule = build_conditions(summaries)

    condition_groups = {g for s in summaries if s.is_core and s.contributes for g in s.groups}
    if manual is not None and manual.conditions is not None:
        rule = replace(rule, conditions=manual.conditions, condition_terms=[manual.conditions],
                       repeated_domains={}, duplicate_condition_groups={}, redundant_with_conditions=[])
        wanted = _identifiers(manual.conditions) - _negated_identifiers(manual.conditions)
        condition_groups = {g for s in summaries if s.all_signature_domains & wanted for g in s.groups}
    if manual is not None and manual.extenders is not None:
        rule = replace(rule, extenders=None if manual.extenders == "-" else manual.extenders, redundant_extenders=[])

    spans: list[ScSpan] = []
    cutoff = None
    if genbank_dir.is_dir():
        memberships = load_groups_tsv(groups_tsv)
        spans = compute_sc_spans(memberships, genbank_dir, condition_groups)
        cutoff = estimate_cutoff_kb(spans)
    if manual is not None and manual.cutoff is not None:
        cutoff = manual.cutoff

    description = product_names.get(rule_id)

    core = [s for s in summaries if s.is_core]
    extra = [s for s in summaries if not s.is_core]
    unresolved_scs = [s.sc_id for s in spans if s.span_kb is None]
    partial_scs = [s for s in spans if s.span_kb is not None and s.unmatched]
    known_domains = {d for s in summaries for d in s.all_signature_domains}
    warnings = build_warnings(
        core, extra, rule, unresolved_scs, partial_scs, description, rule_id, building_blocks_json,
        manual, known_domains, searched_models(segment_dir.parent),
    )

    text = render(
        segment_dir.name, summaries, rule, spans, cutoff, genbank_dir,
        rule_id, description, warnings, manual,
    )
    (segment_dir / "rule_conditions.md").write_text(text)

    rule_block_lines = render_rule_block_with_warnings(
        rule_id, description or "<unknown>", cutoff, rule, warnings, manual,
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

    product_names = load_product_names(args.building_blocks_json) if args.building_blocks_json.exists() else {}
    try:
        manual_rules = load_manual_rules(args.manual_rules) if args.manual_rules.exists() else {}
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    known_ids = {d.name.split("_", 1)[0] for d in group_dirs}
    unknown_ids = sorted(set(manual_rules) - known_ids)
    if unknown_ids and not single_group:
        print(f"warning: {args.manual_rules} has rule(s) with no group directory: {', '.join(unknown_ids)}",
              file=sys.stderr)

    failed = []
    with_warnings: list[str] = []
    with_notes: list[str] = []
    dropped: list[str] = []
    manual_fields: list[str] = []
    for group_dir in group_dirs:
        genbank_dir = args.genbank_dir or (args.genbank_root / group_dir.name)
        try:
            text, warnings = build_group_rule(group_dir, genbank_dir, args.building_blocks_json, product_names,
                                              manual_rules)
        except FileNotFoundError as e:
            print(f"error: {e}", file=sys.stderr)
            failed.append(group_dir.name)
            continue

        rule_id = group_dir.name.split("_", 1)[0]
        manual = manual_rules.get(rule_id)
        if manual is not None and manual.drop:
            dropped.append(rule_id)
        else:
            blocking, advisory = split_advisory(warnings)
            if blocking:
                with_warnings.append(rule_id)
            if advisory:
                with_notes.append(rule_id)
            if manual is not None:
                manual_fields.append(rule_id)

        if single_group:
            print(text)
            print(f"Written to {group_dir / 'rule_conditions.md'} and {group_dir / 'rule.txt'} "
                  f"({len(warnings)} warning(s))", file=sys.stderr)
        else:
            blocking, advisory = split_advisory(warnings)
            bits = ([f"{len(blocking)} warning(s)"] if blocking else []) + ([f"{len(advisory)} note(s)"] if advisory else [])
            status = ", ".join(bits) or "ok"
            if manual is not None:
                status += ", DROPPED (manual)" if manual.drop else ", manual fields"
            print(f"  {group_dir.name}: {status}")

    if not single_group:
        print(f"Done: {len(group_dirs) - len(failed)} group(s) built, {len(failed)} failed.")
        print(f"  {len(with_warnings)} with warnings (skipped by combine_rules.py)"
              + (f": {', '.join(with_warnings)}" if with_warnings else ""))
        print(f"  {len(with_notes)} with advisory notes (still combined)"
              + (f": {', '.join(with_notes)}" if with_notes else ""))
        print(f"  {len(dropped)} dropped in {args.manual_rules.name}" + (f": {', '.join(dropped)}" if dropped else ""))
        print(f"  {len(manual_fields)} with hand-written fields" + (f": {', '.join(manual_fields)}" if manual_fields else ""))
        n_clean = len(group_dirs) - len(failed) - len(with_warnings) - len(dropped)
        print(f"  {n_clean} ready for combine_rules.py")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
