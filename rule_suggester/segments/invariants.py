"""Invariant checks (spec step 7). Must run every time; violations are collected, not silently
swallowed -- callers decide whether to raise.
"""

from __future__ import annotations

from collections import defaultdict

from .grouping import Group
from .io_fasta import Protein
from .segmentation import Segment


def check_invariants(
    proteins: list[Protein],
    segments: list[Segment],
    groups: list[Group],
) -> list[str]:
    violations: list[str] = []

    segs_by_protein: dict[str, list[Segment]] = defaultdict(list)
    for seg in segments:
        segs_by_protein[seg.key].append(seg)

    for p in proteins:
        segs = sorted(segs_by_protein[p.key], key=lambda s: s.start)
        if not segs:
            violations.append(f"protein {p.key}: no segments")
            continue
        if segs[0].start != 0 or segs[-1].end != p.length:
            violations.append(
                f"protein {p.key}: segments span {segs[0].start}-{segs[-1].end}, "
                f"expected 0-{p.length}"
            )
        for a, b in zip(segs, segs[1:]):
            if a.end != b.start:
                violations.append(
                    f"protein {p.key}: segments {a.id} and {b.id} are not contiguous "
                    f"({a.end} != {b.start})"
                )

    segment_ids = {seg.id for seg in segments}
    group_of_segment: dict[str, str] = {}
    for group in groups:
        for seg in group.segments:
            if seg.id in group_of_segment:
                violations.append(
                    f"segment {seg.id} is in groups {group_of_segment[seg.id]} and {group.id}"
                )
            group_of_segment[seg.id] = group.id

    for seg_id in segment_ids:
        if seg_id not in group_of_segment:
            violations.append(f"segment {seg_id} is in no group")

    proteins_with_key = {p.key for p in proteins}
    for group in groups:
        seen_scs: dict[str, str] = {}
        for seg in group.segments:
            if seg.key not in proteins_with_key:
                violations.append(f"group {group.id}: segment {seg.id} has no source protein")
            if seg.sc_id in seen_scs:
                violations.append(
                    f"group {group.id}: two segments from subcluster {seg.sc_id} "
                    f"({seen_scs[seg.sc_id]} and {seg.id})"
                )
            seen_scs[seg.sc_id] = seg.id

    return violations
