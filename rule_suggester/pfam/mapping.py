"""Mapping domain envelopes onto segments, and segments onto group-level domain prevalence
(spec Phase 2, steps 2-3).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, replace

from .clans import Clan
from .params import PfamParams
from .phase1_io import GroupMembership, SegmentRecord
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
        segs_by_protein[(seg.bgc_id, seg.protein_id)].append(seg)

    assignments = []
    for dom in domains:
        candidates = segs_by_protein.get((dom.bgc_id, dom.protein_id), [])
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
    seg_by_coords = {(s.bgc_id, s.protein_id, s.start, s.end): s.segment_id for s in segments}
    group_of_segment = {}
    for g in groups:
        seg_id = seg_by_coords.get((g.bgc_id, g.protein_id, g.start, g.end))
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
    n_bgcs_with_domain: int
    n_bgcs_total: int
    prevalence: float
    is_signature: bool
    any_straddling: bool
    best_i_evalue: float
    is_representative: bool = False


def compute_group_domains(
    assignments: list[DomainAssignment],
    group_of_segment: dict[str, str],
    n_bgcs_by_group: dict[str, int],
    params: PfamParams,
    clan_map: dict[str, Clan] | None = None,
) -> list[GroupDomain]:
    clan_map = clan_map or {}
    bgcs_with_domain: dict[tuple[str, str], set[str]] = defaultdict(set)
    any_straddle: dict[tuple[str, str], bool] = defaultdict(bool)
    best_i_evalue: dict[tuple[str, str], float] = {}
    pfam_names: dict[str, str] = {}

    for a in assignments:
        if a.segment_id is None:
            continue
        group_id = group_of_segment.get(a.segment_id)
        if group_id is None:
            continue
        key = (group_id, a.domain.pfam_acc)
        pfam_names[a.domain.pfam_acc] = a.domain.pfam_name
        if a.straddling:
            any_straddle[key] = True
        else:
            bgcs_with_domain[key].add(a.domain.bgc_id)
            if a.domain.i_evalue < best_i_evalue.get(key, float("inf")):
                best_i_evalue[key] = a.domain.i_evalue

    results = []
    for key in set(bgcs_with_domain) | set(any_straddle):
        group_id, pfam_acc = key
        n_total = n_bgcs_by_group.get(group_id, 0)
        n_with = len(bgcs_with_domain.get(key, set()))
        prevalence = n_with / n_total if n_total else 0.0
        straddle = any_straddle.get(key, False)
        clan = clan_map.get(pfam_acc)
        results.append(
            GroupDomain(
                group_id=group_id,
                pfam_acc=pfam_acc,
                pfam_name=pfam_names.get(pfam_acc, ""),
                clan_acc=clan.acc if clan else "",
                clan_name=clan.name if clan else "",
                n_bgcs_with_domain=n_with,
                n_bgcs_total=n_total,
                prevalence=prevalence,
                is_signature=(prevalence >= params.min_prevalence and not straddle),
                any_straddling=straddle,
                best_i_evalue=best_i_evalue.get(key, float("inf")),
            )
        )
    return _select_representatives(results)


def _select_representatives(rows: list[GroupDomain]) -> list[GroupDomain]:
    """Within each group, nominate one representative Pfam per clan: the signature domain with
    the best (lowest) independent E-value. Signature Pfams with no clan assignment compete with
    nobody, so each is its own representative.
    """
    best_by_clan: dict[tuple[str, str], GroupDomain] = {}
    for g in rows:
        if not g.is_signature:
            continue
        clan_key = g.clan_acc or g.pfam_acc
        key = (g.group_id, clan_key)
        current = best_by_clan.get(key)
        if current is None or g.best_i_evalue < current.best_i_evalue:
            best_by_clan[key] = g

    representative_ids = {(g.group_id, g.pfam_acc) for g in best_by_clan.values()}
    return [replace(g, is_representative=(g.group_id, g.pfam_acc) in representative_ids) for g in rows]
