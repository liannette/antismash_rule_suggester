"""Building antiSMASH-style rule CONDITIONS/EXTENDERS text from representative Pfam domains
(Phase 3, building on Phase 2's per-group representative-domain selection in rule_suggester/pfam).

CONDITIONS: one cds(...) per core group (present in every subcluster of the product group), ANDed --
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

Complementary groups: phase 1 sometimes splits one function into non-core groups that never
co-occur but together cover every subcluster (e.g. a divergent ortholog in one example). If they
share a signature domain, they're merged into one core group (id "G03+G04") whose representatives
are the shared signature domains -- see merge_complementary_groups().

This module only builds the CONDITIONS/EXTENDERS fragments and the per-group detail behind them --
not a complete RULE block, since CATEGORY/CUTOFF/NEIGHBOURHOOD aren't derivable from segment/domain
data and a placeholder would misrepresent them as decided.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass, field, replace
from itertools import combinations
from pathlib import Path

from ..phase1_io import load_groups_tsv


@dataclass(frozen=True)
class _GroupDomains:
    representatives: dict[str, list[str]]  # group_id -> sorted representative names
    signatures: dict[str, set[str]]  # group_id -> signature names
    evalues: dict[str, dict[str, float]]  # group_id -> name -> best i-Evalue
    clans: dict[str, str]  # name -> clan accession (Pfam); absent = no clan


def _load_group_domains(group_domains_tsv: Path) -> _GroupDomains:
    reps: dict[str, list[str]] = defaultdict(list)
    sigs: dict[str, set[str]] = defaultdict(set)
    evalues: dict[str, dict[str, float]] = defaultdict(dict)
    clans: dict[str, str] = {}
    with Path(group_domains_tsv).open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            g, name = row["group_id"], row["pfam_name"]
            if row["is_representative"] == "True":
                reps[g].append(name)
            if row["is_signature"] == "True":
                sigs[g].add(name)
            evalues[g][name] = float(row["best_i_evalue"])
            if row.get("clan_acc"):
                clans[name] = row["clan_acc"]
    for names in reps.values():
        names.sort()
    return _GroupDomains(dict(reps), dict(sigs), dict(evalues), clans)


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
    n_scs: int  # this group's own subcluster membership count
    n_scs_total: int  # total subclusters in the product group
    representative_domains: list[str] = field(default_factory=list)
    signature_domains: frozenset[str] = frozenset()
    domain_evalues: dict[str, float] = field(default_factory=dict)  # best i-Evalue per domain
    member_groups: tuple[str, ...] = ()  # phase-1 groups merged into this one; empty = just group_id
    member_signature_domains: frozenset[str] = frozenset()  # merged group: every member's signatures

    @property
    def all_signature_domains(self) -> frozenset[str]:
        """Every signature domain of the phase-1 groups behind this summary (for a merged group not
        just the shared ones its representatives come from)."""
        return self.member_signature_domains or self.signature_domains

    @property
    def groups(self) -> tuple[str, ...]:
        """The phase-1 group ids this summary stands for (several for a merged group)."""
        return self.member_groups or (self.group_id,)

    @property
    def is_core(self) -> bool:
        return self.n_scs_total > 0 and self.n_scs == self.n_scs_total

    @property
    def contributes(self) -> bool:
        return bool(self.representative_domains)


def _load_presence(presence_tsv: Path) -> tuple[int, dict[str, int], dict[str, set[str]]]:
    """Total subcluster count for the product group (row count), and each group's own subcluster count
    (column sum) -- a group is core iff its own count equals the total. groups.tsv alone can't
    give the total: if no group happens to be core (e.g. a product group whose subclusters share no gene), every group's
    own n_scs would look like a plausible-but-wrong total.
    """
    with Path(presence_tsv).open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        group_ids = next(reader)[1:]
        counts = dict.fromkeys(group_ids, 0)
        scs: dict[str, set[str]] = {g: set() for g in group_ids}
        n_scs_total = 0
        for row in reader:
            n_scs_total += 1
            for group_id, value in zip(group_ids, row[1:]):
                counts[group_id] += int(value)
                if int(value):
                    scs[group_id].add(row[0])
    return n_scs_total, counts, scs


def build_group_summaries(
    groups_tsv: Path, presence_tsv: Path, group_domains_tsv: Path, merge_complementary: bool = True
) -> list[GroupSummary]:
    n_scs_total, group_sc_counts, group_scs = _load_presence(presence_tsv)

    gene_names: dict[str, list[str]] = {}
    order: list[str] = []
    for m in load_groups_tsv(groups_tsv):
        if m.group_id not in gene_names:
            order.append(m.group_id)
            gene_names[m.group_id] = []
        if m.gene_name and m.gene_name not in gene_names[m.group_id]:
            gene_names[m.group_id].append(m.gene_name)

    domains = _load_group_domains(group_domains_tsv)

    summaries = [
        GroupSummary(
            group_id=g,
            gene_names=gene_names[g],
            n_scs=group_sc_counts.get(g, 0),
            n_scs_total=n_scs_total,
            representative_domains=domains.representatives.get(g, []),
            signature_domains=frozenset(domains.signatures.get(g, set())),
            domain_evalues=domains.evalues.get(g, {}),
        )
        for g in order
    ]
    if not merge_complementary:
        return summaries
    all_scs = set().union(*group_scs.values()) if group_scs else set()
    return merge_complementary_groups(summaries, group_scs, all_scs, domains.clans)


def merge_complementary_groups(
    summaries: list[GroupSummary], group_scs: dict[str, set[str]], all_scs: set[str],
    clans: dict[str, str] | None = None,
) -> list[GroupSummary]:
    """Merge non-core groups that never share a subcluster, together cover all of them, and share
    at least one signature domain, into one core group. Its representatives are the shared signature
    domains that are representative in at least one member; if there are none, the shared signature
    domains with one per Pfam clan (the best worst-case i-Evalue across members), as step 2 does.
    Candidate sets are tried smallest first; a group joins at most one merge.
    """
    clans = clans or {}
    noncore = [s for s in summaries if not s.is_core and s.signature_domains]
    by_id = {s.group_id: s for s in summaries}
    merges: list[tuple[str, ...]] = []
    used: set[str] = set()
    for size in range(2, len(noncore) + 1):
        for combo in combinations(noncore, size):
            ids = tuple(s.group_id for s in combo)
            if used & set(ids):
                continue
            covered = [group_scs.get(i, set()) for i in ids]
            if sum(len(c) for c in covered) != len(set().union(*covered)) or set().union(*covered) != all_scs:
                continue
            if frozenset.intersection(*(s.signature_domains for s in combo)):
                merges.append(ids)
                used |= set(ids)
    if not merges:
        return summaries

    merged_by_first: dict[str, GroupSummary] = {}
    for ids in merges:
        members = [by_id[i] for i in ids]
        shared = frozenset.intersection(*(m.signature_domains for m in members))
        reps = sorted(shared & {d for m in members for d in m.representative_domains})
        if not reps:
            best_by_clan: dict[str, tuple[float, str]] = {}
            for d in sorted(shared):
                worst = max(m.domain_evalues.get(d, 1.0) for m in members)
                key = clans.get(d, d)
                if key not in best_by_clan or worst < best_by_clan[key][0]:
                    best_by_clan[key] = (worst, d)
            reps = sorted(d for _, d in best_by_clan.values())
        merged_by_first[ids[0]] = GroupSummary(
            group_id="+".join(ids),
            gene_names=[g for m in members for g in m.gene_names],
            n_scs=len(all_scs),
            n_scs_total=members[0].n_scs_total,
            representative_domains=reps,
            signature_domains=shared,
            domain_evalues={d: max(m.domain_evalues.get(d, 1.0) for m in members) for d in shared},
            member_groups=ids,
            member_signature_domains=frozenset().union(*(m.signature_domains for m in members)),
        )
    out = []
    for s in summaries:
        if s.group_id in merged_by_first:
            out.append(merged_by_first[s.group_id])
        elif s.group_id not in used:
            out.append(s)
    return out


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
