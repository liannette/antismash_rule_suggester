"""Building antiSMASH-style rule CONDITIONS/EXTENDERS text from representative Pfam domains
(Phase 3, building on Phase 2's per-group representative-domain selection in rule_suggester/pfam).

CONDITIONS: one cds(...) per core group (present in every BGC of the product group), ANDed --
a core group's representative domains must co-occur on the same CDS, and every core group must
be hit for the rule to fire. A core group with no signature domain at all (nothing reached 100%
prevalence in Step 2) has nothing to contribute and is left out of CONDITIONS -- that weakens the
rule silently unless it's surfaced, so callers should report it (see GroupSummary.contributes).

When several core groups are paralogous and each nominates the very same Pfam as its sole
representative (e.g. three crotonase-fold genes all represented by the same fold's dominant
profile), a plain 'and' of that id N times would just be deduped to one occurrence -- antiSMASH
checks profile presence, not hit count, so a bare id can't require N separate genes. repeatable(N,
[id]) does: it counts hits of that profile across nearby CDSs rather than requiring N distinct
profiles, so it correctly encodes "this domain must be hit on N different genes". Only single-
domain (singleton) group terms are eligible -- repeatable()'s option list has no way to require a
specific combination of *several* domains to co-occur per repeat, so multi-domain group terms
keep their own cds() term instead. repeatable() is also only legal directly in CONDITIONS (the
antiSMASH grammar disallows it inside cds() and inside EXTENDERS), which is exactly where this
applies.

EXTENDERS: representative domains of non-core (accessory) groups, ORed together inside a single
cds(...). antiSMASH's EXTENDERS field only accepts one flat, unparenthesised cds()-scoped
condition -- no nested cds()/parens -- which is exactly what an OR-of-AND-chains gives via the
DSL's 'and' > 'or' precedence: each group's own AND-chain stays bound together, and the groups
are alternatives to each other.

Rule identifiers must be the Pfam NAME (e.g. "GFO_IDH_MocA"), not the PFxxxxx accession -- that's
what antiSMASH's HMM library indexes profiles by.

This module only builds the CONDITIONS/EXTENDERS fragments and the per-group detail behind them --
not a complete RULE block, since CATEGORY/CUTOFF/NEIGHBOURHOOD aren't derivable from segment/domain
data and a placeholder would misrepresent them as decided.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from ..phase1_io import load_groups_tsv


def load_representative_domains(group_domains_tsv: Path) -> dict[str, list[str]]:
    """group_id -> sorted list of representative Pfam names."""
    reps: dict[str, list[str]] = {}
    with Path(group_domains_tsv).open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["is_representative"] == "True":
                reps.setdefault(row["group_id"], []).append(row["pfam_name"])
    for names in reps.values():
        names.sort()
    return reps


@dataclass(frozen=True)
class GroupSummary:
    group_id: str
    gene_names: list[str]
    n_bgcs: int  # this group's own BGC membership count
    n_bgcs_total: int  # total BGCs in the product group
    representative_domains: list[str] = field(default_factory=list)

    @property
    def is_core(self) -> bool:
        return self.n_bgcs_total > 0 and self.n_bgcs == self.n_bgcs_total

    @property
    def contributes(self) -> bool:
        return bool(self.representative_domains)


def _load_presence(presence_tsv: Path) -> tuple[int, dict[str, int]]:
    """Total BGC count for the product group (row count), and each group's own BGC count
    (column sum) -- a group is core iff its own count equals the total. groups.tsv alone can't
    give the total: if no group happens to be core (e.g. SCG0045_ethylmalonyl-CoA), every group's
    own n_bgcs would look like a plausible-but-wrong total.
    """
    with Path(presence_tsv).open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        group_ids = next(reader)[1:]
        counts = dict.fromkeys(group_ids, 0)
        n_bgcs_total = 0
        for row in reader:
            n_bgcs_total += 1
            for group_id, value in zip(group_ids, row[1:]):
                counts[group_id] += int(value)
    return n_bgcs_total, counts


def build_group_summaries(groups_tsv: Path, presence_tsv: Path, group_domains_tsv: Path) -> list[GroupSummary]:
    n_bgcs_total, group_bgc_counts = _load_presence(presence_tsv)

    gene_names: dict[str, list[str]] = {}
    order: list[str] = []
    for m in load_groups_tsv(groups_tsv):
        if m.group_id not in gene_names:
            order.append(m.group_id)
            gene_names[m.group_id] = []
        if m.gene_name and m.gene_name not in gene_names[m.group_id]:
            gene_names[m.group_id].append(m.gene_name)

    reps = load_representative_domains(group_domains_tsv)

    return [
        GroupSummary(
            group_id=g,
            gene_names=gene_names[g],
            n_bgcs=group_bgc_counts.get(g, 0),
            n_bgcs_total=n_bgcs_total,
            representative_domains=reps.get(g, []),
        )
        for g in order
    ]


def _group_condition(pfam_names: list[str]) -> str:
    """One group's CONDITIONS term: a bare id if it has a single representative domain, else a
    cds() AND-chain requiring them all on the same protein.
    """
    if len(pfam_names) == 1:
        return pfam_names[0]
    return f"cds({' and '.join(pfam_names)})"


@dataclass(frozen=True)
class RuleConditions:
    conditions: str | None  # None if no core group has a representative domain
    extenders: str | None  # None if no non-core group has a representative domain
    condition_terms: list[str] = field(default_factory=list)  # the same terms, unjoined -- for multi-line rendering
    repeated_domains: dict[str, list[str]] = field(default_factory=dict)  # domain -> group_ids folded into repeatable()
    duplicate_condition_groups: dict[str, list[str]] = field(default_factory=dict)  # multi-domain term -> group_ids
    redundant_extenders: list[tuple[str, str]] = field(default_factory=list)  # (dominated_group, dominating_group)
    redundant_with_conditions: list[str] = field(default_factory=list)  # non-core group_ids dropped from EXTENDERS


def _core_terms(core: list[GroupSummary]) -> tuple[list[str], dict[str, list[str]], dict[str, list[str]]]:
    """Core groups whose sole representative domain is shared with other core groups become one
    repeatable(N, [domain]) term instead of N redundant copies of the same bare id; every other
    group keeps its own term, deduped if identical (e.g. two different core groups both nominating
    the exact same multi-domain combination -- reported via duplicate_condition_groups, since that
    can't be folded into repeatable() the way a single shared domain can).
    """
    singleton_groups: dict[str, list[str]] = defaultdict(list)
    other_group_ids_by_term: dict[str, list[str]] = defaultdict(list)
    for s in core:
        if len(s.representative_domains) == 1:
            singleton_groups[s.representative_domains[0]].append(s.group_id)
        else:
            other_group_ids_by_term[_group_condition(s.representative_domains)].append(s.group_id)

    repeated_domains: dict[str, list[str]] = {}
    singleton_terms: list[str] = []
    for domain, group_ids in singleton_groups.items():
        if len(group_ids) >= 2:
            singleton_terms.append(f"repeatable({len(group_ids)}, [{domain}])")
            repeated_domains[domain] = group_ids
        else:
            singleton_terms.append(domain)

    duplicate_condition_groups = {
        term: group_ids for term, group_ids in other_group_ids_by_term.items() if len(group_ids) >= 2
    }

    terms = list(dict.fromkeys(singleton_terms + list(other_group_ids_by_term)))
    return terms, repeated_domains, duplicate_condition_groups


def _redundant_extenders(extra: list[GroupSummary]) -> list[tuple[str, str]]:
    """EXTENDERS ORs each non-core group's AND-chain together. If one group's domain set is a
    strict superset of another's, the superset group's clause is logically dead: whenever it's
    satisfied, the subset group's (weaker) clause is already satisfied too, so the OR never needs
    it. This doesn't change what gets written to EXTENDERS (removing it could hide a genuine
    multi-domain signature that's independently informative) -- it's reported so it can be
    reviewed by hand.
    """
    redundant = []
    for s in extra:
        domains = set(s.representative_domains)
        for other in extra:
            if other.group_id != s.group_id and domains > set(other.representative_domains):
                redundant.append((s.group_id, other.group_id))
                break
    return redundant


def build_conditions(summaries: list[GroupSummary]) -> RuleConditions:
    """Two different groups can nominate the same Pfam (e.g. two aminotransferase paralogs each
    represent their clan by the same profile) -- see module docstring for why that becomes
    repeatable() rather than a deduped-away or meaninglessly repeated bare id.
    """
    core = [s for s in summaries if s.is_core and s.contributes]
    extra = [s for s in summaries if not s.is_core and s.contributes]

    core_terms, repeated_domains, duplicate_condition_groups = _core_terms(core)
    conditions = " and ".join(core_terms) if core_terms else None

    # A non-core group whose every representative domain is already required by CONDITIONS (e.g.
    # an accessory paralog nominating the same Pfam a core group already anchors on) contributes
    # nothing to EXTENDERS: that OR-term is guaranteed satisfied whenever CONDITIONS is, making the
    # whole EXTENDERS clause a no-op. Drop it and report it instead of silently building a rule
    # that looks like it extends the cluster but never can.
    core_domains = {d for s in core for d in s.representative_domains}
    redundant_with_conditions = [s.group_id for s in extra if set(s.representative_domains) <= core_domains]
    extra = [s for s in extra if not set(s.representative_domains) <= core_domains]

    extenders = None
    if extra:
        extra_chains = dict.fromkeys(" and ".join(s.representative_domains) for s in extra)
        chain = " or ".join(extra_chains)
        # antiSMASH's EXTENDERS grammar accepts a bare identifier or a cds(...) of 2+ identifiers,
        # never cds() around a single one -- so only skip the wrapper when every non-core group
        # reduces (after dedup) to the exact same single domain.
        single_domain = len(extra_chains) == 1 and " and " not in next(iter(extra_chains))
        extenders = chain if single_domain else f"cds({chain})"

    return RuleConditions(
        conditions=conditions,
        extenders=extenders,
        condition_terms=core_terms,
        repeated_domains=repeated_domains,
        duplicate_condition_groups=duplicate_condition_groups,
        redundant_extenders=_redundant_extenders(extra),
        redundant_with_conditions=redundant_with_conditions,
    )
