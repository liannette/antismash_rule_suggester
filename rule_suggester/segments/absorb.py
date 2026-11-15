"""Collapsing breakpoints that don't discriminate between subclusters.

Two adjacent same-protein segments landing in different groups only matters if their subcluster coverage
actually differs -- i.e. some subcluster has one piece but not the other. If the observed same-protein
adjacency between two groups covers *every* member of at least one of them, the breakpoint between
them carries no such distinguishing information (every subcluster on that side already has both pieces
back to back), so the groups are merged: boundary extended, breakpoint dropped. This covers both:
  - a small/partial-coverage segment folded into a much bigger, better-supported neighbour
    (adjacency covers all of the *smaller* group), and
  - two equally-supported segments that are simply one domain split in two by HSP fragmentation,
    with identical subcluster coverage (adjacency covers all of *both* groups).

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
        cluster_scs: dict[str, set[str]] = defaultdict(set)
        for seg in segments_by_id.values():
            cluster_scs[cluster_of[seg.id]].add(seg.sc_id)

        by_protein: dict[str, list[Segment]] = defaultdict(list)
        for seg in segments_by_id.values():
            by_protein[seg.key].append(seg)
        for key in by_protein:
            by_protein[key].sort(key=lambda s: s.start)

        # (root_a, root_b) -> subclusters where a's protein has segment a directly followed by segment b
        adjacency: dict[tuple[str, str], set[str]] = defaultdict(set)
        for segs in by_protein.values():
            for a, b in zip(segs, segs[1:]):
                if a.end != b.start:
                    continue  # not contiguous -- leave any gap/overlap alone
                ra, rb = cluster_of[a.id], cluster_of[b.id]
                if ra != rb:
                    adjacency[(ra, rb)].add(a.sc_id)

        target = None
        for pair in sorted(adjacency):
            ra, rb = pair
            support = adjacency[pair]
            if support == cluster_scs[ra] or support == cluster_scs[rb]:
                target = pair
                break
        if target is None:
            break

        ra, rb = target
        # The merged groups continue as the one whose membership the adjacency does NOT fully
        # cover (the bigger, better-supported group); the fully covered one is absorbed into it.
        # Keeping `ra` unconditionally would instead pull a homolog out of its real group into
        # the singleton of its own N-terminal fragment.
        support = adjacency[target]
        survivor, absorbed = (rb, ra) if support == cluster_scs[ra] and support != cluster_scs[rb] else (ra, rb)
        for segs in by_protein.values():
            for a, b in zip(segs, segs[1:]):
                if cluster_of.get(a.id) == ra and cluster_of.get(b.id) == rb:
                    protein = proteins_by_key[a.key]
                    pid = protein.protein_id or protein.local_key
                    new_seg = Segment(
                        id=f"{protein.sc_id}|{pid}|{a.start + 1}-{b.end}",
                        key=a.key,
                        sc_id=a.sc_id,
                        protein_id=a.protein_id,
                        gene_name=a.gene_name,
                        locus_tag=a.locus_tag,
                        start=a.start,
                        end=b.end,
                    )
                    del segments_by_id[a.id], segments_by_id[b.id]
                    del cluster_of[a.id], cluster_of[b.id]
                    segments_by_id[new_seg.id] = new_seg
                    cluster_of[new_seg.id] = survivor
                    merges.append(SegmentMerge(a.key, f"{a.id} + {b.id}", new_seg.id))
                    break  # a protein has at most one (ra,rb)-matching pair
        # any member of the absorbed group not part of a merged pair joins the survivor too
        for seg_id, root in cluster_of.items():
            if root == absorbed:
                cluster_of[seg_id] = survivor

    return list(segments_by_id.values()), cluster_of, merges
