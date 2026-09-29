"""Assigning final group IDs to clusters of segments.

Order: core groups (present in all BGCs) first, then by descending n_bgcs, then by mean position
of the first member gene within its BGC's FASTA order; ties broken by the lowest member segment id
for determinism. IDs are G01, G02, ... (zero-padded to fit the group count).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from statistics import mean

from .io_fasta import Protein
from .segmentation import Segment


@dataclass(frozen=True)
class Group:
    id: str
    segments: list[Segment]

    @property
    def n_bgcs(self) -> int:
        return len({seg.bgc_id for seg in self.segments})


def assign_groups(
    segments: list[Segment],
    cluster_of: dict[str, str],
    proteins: list[Protein],
    total_bgcs: int,
) -> list[Group]:
    order = {p.key: p.order for p in proteins}

    clusters: dict[str, list[Segment]] = defaultdict(list)
    for seg in segments:
        clusters[cluster_of[seg.id]].append(seg)

    def sort_key(members: list[Segment]):
        n_bgcs = len({seg.bgc_id for seg in members})
        is_core = n_bgcs == total_bgcs
        mean_pos = mean(order[seg.key] for seg in members)
        min_id = min(seg.id for seg in members)
        return (not is_core, -n_bgcs, mean_pos, min_id)

    ordered_clusters = sorted(clusters.values(), key=sort_key)

    width = max(2, len(str(len(ordered_clusters))))
    groups = []
    for i, members in enumerate(ordered_clusters, start=1):
        members_sorted = sorted(members, key=lambda seg: (seg.bgc_id, seg.start))
        groups.append(Group(id=f"G{i:0{width}d}", segments=members_sorted))
    return groups
