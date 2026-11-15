"""Mapping domain envelopes onto segments, and segments onto group-level domain prevalence
(spec Phase 2, steps 2-3).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, replace

from .clans import Clan
from .params import PfamParams
from ..phase1_io import GroupMembership, SegmentRecord
from .scan import Domain


@dataclass(frozen=True)
class DomainAssignment:
    domain: Domain
    segment_id: str | None  # None if the domain's protein has no matching segment in this group
    frac_domain_in_segment: float  # share of the domain's own length landing in that segment
    frac_segment_covered: float  # share of the segment's length the domain covers
    straddling: bool  # no single segment reaches straddle_threshold of the domain


def map_domains_to_segments(
    domains: list[Domain], segments: list[SegmentRecord], params: PfamParams
) -> list[DomainAssignment]:
    segs_by_protein: dict[tuple[str, str], list[SegmentRecord]] = defaultdict(list)
    for seg in segments:
        segs_by_protein[(seg.sc_id, seg.protein_id)].append(seg)

    assignments = []
    for dom in domains:
        candidates = segs_by_protein.get((dom.sc_id, dom.protein_id), [])
        domain_len = dom.env_to - dom.env_from

        best_seg, best_frac_domain, best_frac_segment = None, 0.0, 0.0
        for seg in candidates:
            ov_start, ov_end = max(dom.env_from, seg.start), min(dom.env_to, seg.end)
            if ov_end <= ov_start:
                continue
            overlap = ov_end - ov_start
            frac_domain = overlap / domain_len if domain_len > 0 else 0.0
            if frac_domain > best_frac_domain:
                best_seg, best_frac_domain = seg, frac_domain
                best_frac_segment = overlap / seg.length if seg.length > 0 else 0.0

        assignments.append(
            DomainAssignment(
                domain=dom,
                segment_id=best_seg.segment_id if best_seg else None,
                frac_domain_in_segment=best_frac_domain,
                frac_segment_covered=best_frac_segment,
                straddling=best_frac_domain < params.straddle_threshold,
            )
        )
    return assignments


def build_group_of_segment(
    segments: list[SegmentRecord], groups: list[GroupMembership]
) -> dict[str, str]:
    """groups.tsv doesn't carry segment_id, so join it back to segments.tsv by coordinates."""
    seg_by_coords = {(s.sc_id, s.protein_id, s.start, s.end): s.segment_id for s in segments}
    group_of_segment = {}
    for g in groups:
        seg_id = seg_by_coords.get((g.sc_id, g.protein_id, g.start, g.end))
        if seg_id is not None:
            group_of_segment[seg_id] = g.group_id
    return group_of_segment


@dataclass(frozen=True)
class GroupDomain:
    group_id: str
    pfam_acc: str
    pfam_name: str
    clan_acc: str
    clan_name: str
    n_scs_with_domain: int
    n_scs_total: int
    prevalence: float
    is_signature: bool
    any_straddling: bool
    best_i_evalue: float
    is_representative: bool = False
    source: str = "Pfam"
    family_type: str = ""
    description: str = ""
    replaced_by: str = ""  # signature not representative because this PGAP representative covers it


Envelopes = dict[tuple[str, str], dict[str, list[tuple[int, int]]]]  # (group, acc) -> sc -> envelopes


def compute_group_domains(
    assignments: list[DomainAssignment],
    group_of_segment: dict[str, str],
    n_scs_by_group: dict[str, int],
    params: PfamParams,
    clan_map: dict[str, Clan] | None = None,
) -> list[GroupDomain]:
    clan_map = clan_map or {}
    scs_with_domain: dict[tuple[str, str], set[str]] = defaultdict(set)
    any_straddle: dict[tuple[str, str], bool] = defaultdict(bool)
    best_i_evalue: dict[tuple[str, str], float] = {}
    first_domain: dict[str, Domain] = {}
    envelopes: Envelopes = defaultdict(lambda: defaultdict(list))

    for a in assignments:
        if a.segment_id is None:
            continue
        group_id = group_of_segment.get(a.segment_id)
        if group_id is None:
            continue
        key = (group_id, a.domain.pfam_acc)
        first_domain.setdefault(a.domain.pfam_acc, a.domain)
        if a.straddling:
            any_straddle[key] = True
        else:
            scs_with_domain[key].add(a.domain.sc_id)
            envelopes[key][a.domain.sc_id].append((a.domain.env_from, a.domain.env_to))
            if a.domain.i_evalue < best_i_evalue.get(key, float("inf")):
                best_i_evalue[key] = a.domain.i_evalue

    results = []
    for key in set(scs_with_domain) | set(any_straddle):
        group_id, pfam_acc = key
        n_total = n_scs_by_group.get(group_id, 0)
        n_with = len(scs_with_domain.get(key, set()))
        prevalence = n_with / n_total if n_total else 0.0
        straddle = any_straddle.get(key, False)
        clan = clan_map.get(pfam_acc)
        d = first_domain[pfam_acc]
        results.append(
            GroupDomain(
                group_id=group_id,
                pfam_acc=pfam_acc,
                pfam_name=d.pfam_name,
                clan_acc=clan.acc if clan else "",
                clan_name=clan.name if clan else "",
                n_scs_with_domain=n_with,
                n_scs_total=n_total,
                prevalence=prevalence,
                is_signature=(prevalence >= params.min_prevalence and not straddle),
                any_straddling=straddle,
                best_i_evalue=best_i_evalue.get(key, float("inf")),
                source=d.source,
                family_type=d.family_type,
                description=d.description,
            )
        )
    return _select_representatives(results, envelopes, params.replace_overlap, params.representative_pgap_sources)


# PGAP family types, most specific first: equivalogs name one function, subfamilies/domains less so
_FAMILY_RANK = {
    "exception": 0, "equivalog": 0, "PfamEq": 1, "hypoth_equivalog": 1, "equivalog_domain": 1, "PfamAutoEq": 1,
    "subfamily": 2, "subfamily_domain": 3, "domain": 3, "repeat": 4,
}


def _covers(envelopes: Envelopes, group_id: str, covering: str, covered: str, min_share: float) -> bool:
    """True if model `covering` overlaps >= min_share of model `covered`'s envelope in a majority of
    the subclusters where the group has both."""
    a, b = envelopes.get((group_id, covering), {}), envelopes.get((group_id, covered), {})
    shared = set(a) & set(b)
    if not shared:
        return False
    n_covered = 0
    for sc in shared:
        if any(
            (min(c_to, d_to) - max(c_from, d_from)) >= min_share * (d_to - d_from)
            for c_from, c_to in a[sc] for d_from, d_to in b[sc]
        ):
            n_covered += 1
    return 2 * n_covered > len(shared)


def _select_representatives(
    rows: list[GroupDomain], envelopes: Envelopes, min_share: float, pgap_sources: tuple[str, ...]
) -> list[GroupDomain]:
    """Within each group, among signature models:
    1. PGAP models from pgap_sources (TIGRFAM, NCBIFAM; not PRK): most specific family type first
       (then best i-Evalue); one overlapping a model already chosen is left out (replaced_by it).
       PGAP models from other sources are never representative and replace nothing.
    2. Pfam: a domain covered by a chosen PGAP representative is left out (replaced_by it).
    3. Remaining Pfam: one representative per clan, the best (lowest) independent E-value; Pfams
       with no clan compete with nobody, so each is its own representative.
    """
    by_group: dict[str, list[GroupDomain]] = defaultdict(list)
    for g in rows:
        by_group[g.group_id].append(g)

    out = []
    for group_id, group_rows in by_group.items():
        signature = [g for g in group_rows if g.is_signature]
        replaced: dict[str, str] = {}
        representatives: set[str] = set()

        chosen_pgap: list[GroupDomain] = []
        for g in sorted((g for g in signature if g.source in pgap_sources),
                        key=lambda g: (_FAMILY_RANK.get(g.family_type, 5), g.best_i_evalue, g.pfam_acc)):
            winner = next((c for c in chosen_pgap
                           if _covers(envelopes, group_id, c.pfam_acc, g.pfam_acc, min_share)
                           or _covers(envelopes, group_id, g.pfam_acc, c.pfam_acc, min_share)), None)
            if winner:
                replaced[g.pfam_acc] = winner.pfam_name
            else:
                chosen_pgap.append(g)
                representatives.add(g.pfam_acc)

        best_by_clan: dict[str, GroupDomain] = {}
        for g in signature:
            if g.source != "Pfam":
                continue
            winner = next((c for c in chosen_pgap if _covers(envelopes, group_id, c.pfam_acc, g.pfam_acc, min_share)), None)
            if winner:
                replaced[g.pfam_acc] = winner.pfam_name
                continue
            clan_key = g.clan_acc or g.pfam_acc
            current = best_by_clan.get(clan_key)
            if current is None or g.best_i_evalue < current.best_i_evalue:
                best_by_clan[clan_key] = g
        representatives |= {g.pfam_acc for g in best_by_clan.values()}

        out += [replace(g, is_representative=g.pfam_acc in representatives, replaced_by=replaced.get(g.pfam_acc, ""))
                for g in group_rows]
    return out
