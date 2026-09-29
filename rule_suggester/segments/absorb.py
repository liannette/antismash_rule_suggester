"""Collapsing breakpoints that don't discriminate between BGCs.

Two adjacent same-protein segments landing in different groups only matters if their BGC coverage
actually differs -- i.e. some BGC has one piece but not the other. If the observed same-protein
adjacency between two groups covers *every* member of at least one of them, the breakpoint between
them carries no such distinguishing information (every BGC on that side already has both pieces
back to back), so the groups are merged: boundary extended, breakpoint dropped. This covers both:
  - a small/partial-coverage segment folded into a much bigger, better-supported neighbour
    (adjacency covers all of the *smaller* group), and
  - two equally-supported segments that are simply one domain split in two by HSP fragmentation,
    with identical BGC coverage (adjacency covers all of *both* groups).

A genuine gene fusion is the opposite case: the two segments are only adjacent in the one fused
lineage, a small fraction of either group's full membership, so that adjacency covers neither
group fully and the groups are left alone -- per the spec's requirement that a fused gene
contributes one segment to each group its parts belong to.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from .io_fasta import Protein
from .segmentation import Segment


@dataclass(frozen=True)
class SegmentMerge:
    protein_key: str
    from_ids: str  # the two segment ids that were merged, for the report
    into_id: str


def merge_uninformative_breakpoints(
    proteins: list[Protein],
    segments: list[Segment],
    cluster_of: dict[str, str],
) -> tuple[list[Segment], dict[str, str], list[SegmentMerge]]:
    segments_by_id = {seg.id: seg for seg in segments}
    cluster_of = dict(cluster_of)
    proteins_by_key = {p.key: p for p in proteins}
    merges: list[SegmentMerge] = []

    while True:
        cluster_bgcs: dict[str, set[str]] = defaultdict(set)
        for seg in segments_by_id.values():
            cluster_bgcs[cluster_of[seg.id]].add(seg.bgc_id)

        by_protein: dict[str, list[Segment]] = defaultdict(list)
        for seg in segments_by_id.values():
            by_protein[seg.key].append(seg)
        for key in by_protein:
            by_protein[key].sort(key=lambda s: s.start)

        # (root_a, root_b) -> BGCs where a's protein has segment a directly followed by segment b
        adjacency: dict[tuple[str, str], set[str]] = defaultdict(set)
        for segs in by_protein.values():
            for a, b in zip(segs, segs[1:]):
                if a.end != b.start:
                    continue  # not contiguous -- leave any gap/overlap alone
                ra, rb = cluster_of[a.id], cluster_of[b.id]
                if ra != rb:
                    adjacency[(ra, rb)].add(a.bgc_id)

        target = None
        for pair in sorted(adjacency):
            ra, rb = pair
            support = adjacency[pair]
            if support == cluster_bgcs[ra] or support == cluster_bgcs[rb]:
                target = pair
                break
        if target is None:
            break

        ra, rb = target
        for segs in by_protein.values():
            for a, b in zip(segs, segs[1:]):
                if cluster_of.get(a.id) == ra and cluster_of.get(b.id) == rb:
                    protein = proteins_by_key[a.key]
                    pid = protein.protein_id or protein.local_key
                    new_seg = Segment(
                        id=f"{protein.bgc_id}|{pid}|{a.start + 1}-{b.end}",
                        key=a.key,
                        bgc_id=a.bgc_id,
                        protein_id=a.protein_id,
                        gene_name=a.gene_name,
                        locus_tag=a.locus_tag,
                        start=a.start,
                        end=b.end,
                    )
                    del segments_by_id[a.id], segments_by_id[b.id]
                    del cluster_of[a.id], cluster_of[b.id]
                    segments_by_id[new_seg.id] = new_seg
                    cluster_of[new_seg.id] = ra
                    merges.append(SegmentMerge(a.key, f"{a.id} + {b.id}", new_seg.id))
                    break  # a protein has at most one (ra,rb)-matching pair

    return list(segments_by_id.values()), cluster_of, merges
