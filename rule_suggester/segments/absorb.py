"""Post-clustering absorption of minor same-protein segments into a bigger neighbour.

A protein that was over-split by breakpoint noise can end up with one piece in a small or
non-core group while an adjacent piece (same protein) lands in a much better-supported group. That
minor piece adds no coverage the bigger neighbouring group doesn't already have for this protein's
BGC, so rather than let it stand as its own group, fold it into the bigger neighbouring segment --
extending that segment's boundary and dropping the breakpoint between them.

Only same-protein neighbours are ever merged, so this can't mix two different proteins into one
segment. A merge only fires when one side's group has strictly more BGCs than the other's; a tie
merges into the earlier (left) segment for determinism. Segments not adjacent to a bigger
same-protein sibling -- e.g. a real singleton domain with no homolog anywhere -- are untouched.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from .io_fasta import Protein
from .segmentation import Segment


@dataclass(frozen=True)
class Absorption:
    protein_key: str
    absorbed_id: str  # the smaller segment's id, which no longer exists after this
    into_id: str  # the resulting (extended) segment's id


def _cluster_bgc_counts(segments_by_id: dict[str, Segment], cluster_of: dict[str, str]) -> dict[str, int]:
    bgcs: dict[str, set[str]] = defaultdict(set)
    for seg in segments_by_id.values():
        bgcs[cluster_of[seg.id]].add(seg.bgc_id)
    return {root: len(s) for root, s in bgcs.items()}


def absorb_minor_segments(
    proteins: list[Protein],
    segments: list[Segment],
    cluster_of: dict[str, str],
) -> tuple[list[Segment], dict[str, str], list[Absorption]]:
    segments_by_id = {seg.id: seg for seg in segments}
    cluster_of = dict(cluster_of)
    proteins_by_key = {p.key: p for p in proteins}
    absorptions: list[Absorption] = []

    by_protein: dict[str, list[str]] = defaultdict(list)
    for seg in segments:
        by_protein[seg.key].append(seg.id)

    for key, seg_ids in by_protein.items():
        protein = proteins_by_key[key]
        pid = protein.protein_id or protein.local_key

        while len(seg_ids) > 1:
            seg_ids = sorted(seg_ids, key=lambda sid: segments_by_id[sid].start)
            sizes = _cluster_bgc_counts(segments_by_id, cluster_of)

            merge_pair = None
            for i in range(len(seg_ids) - 1):
                a, b = segments_by_id[seg_ids[i]], segments_by_id[seg_ids[i + 1]]
                size_a, size_b = sizes[cluster_of[a.id]], sizes[cluster_of[b.id]]
                if size_a != size_b:
                    merge_pair = (a, b) if size_a >= size_b else (b, a)  # tie favours the left segment
                    break
            if merge_pair is None:
                break

            bigger, smaller = merge_pair
            new_start, new_end = min(bigger.start, smaller.start), max(bigger.end, smaller.end)
            new_seg = Segment(
                id=f"{protein.bgc_id}|{pid}|{new_start + 1}-{new_end}",
                key=key,
                bgc_id=protein.bgc_id,
                protein_id=protein.protein_id,
                gene_name=protein.gene_name,
                locus_tag=protein.locus_tag,
                start=new_start,
                end=new_end,
            )

            root = cluster_of[bigger.id]
            del segments_by_id[bigger.id], segments_by_id[smaller.id]
            del cluster_of[bigger.id], cluster_of[smaller.id]
            segments_by_id[new_seg.id] = new_seg
            cluster_of[new_seg.id] = root
            absorptions.append(Absorption(key, smaller.id, new_seg.id))

            seg_ids = [sid for sid in seg_ids if sid not in (bigger.id, smaller.id)] + [new_seg.id]

    return list(segments_by_id.values()), cluster_of, absorptions
