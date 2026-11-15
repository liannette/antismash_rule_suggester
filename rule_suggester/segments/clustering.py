"""Constrained average-linkage clustering of segments into groups (spec step 6).

No synteny bonus in phase 1 -- a synteny-aware tie-breaker on the average-linkage score `L` (e.g.
boosting L for clusters whose members keep the same relative gene order across subclusters) would hook in
right where L is computed in _average_linkage, before the MIN_LINK comparison.
"""

from __future__ import annotations

import heapq

from .params import Params
from .segmentation import Segment


def _average_linkage(
    members_a: set[str], members_b: set[str], pair_weight: dict[tuple[str, str], float]
) -> float:
    total = 0.0
    for a in members_a:
        for b in members_b:
            pair = (a, b) if a < b else (b, a)
            total += pair_weight.get(pair, 0.0)
    return total / (len(members_a) * len(members_b))


def cluster_segments(
    segments: list[Segment], edges: list[tuple[str, str, float]], params: Params
) -> dict[str, str]:
    """Returns {segment_id: cluster_root_id}."""
    parent = {seg.id: seg.id for seg in segments}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    cluster_members: dict[str, set[str]] = {seg.id: {seg.id} for seg in segments}
    cluster_scs: dict[str, set[str]] = {seg.id: {seg.sc_id} for seg in segments}

    pair_weight: dict[tuple[str, str], float] = {}
    for a, b, w in edges:
        pair = (a, b) if a < b else (b, a)
        if w > pair_weight.get(pair, 0.0):
            pair_weight[pair] = w

    heap = [(-w, a, b) for (a, b), w in pair_weight.items()]
    heapq.heapify(heap)

    while heap:
        neg_w, a, b = heapq.heappop(heap)
        w = -neg_w
        ra, rb = find(a), find(b)
        if ra == rb:
            continue
        if cluster_scs[ra] & cluster_scs[rb]:
            continue

        link = _average_linkage(cluster_members[ra], cluster_members[rb], pair_weight)
        if link < params.min_link:
            continue
        if link < w - 1e-9:
            heapq.heappush(heap, (-link, a, b))  # stale score; retry once it resurfaces
            continue

        lo, hi = (ra, rb) if ra < rb else (rb, ra)  # deterministic: smaller id survives as root
        parent[hi] = lo
        cluster_members[lo] |= cluster_members[hi]
        cluster_scs[lo] |= cluster_scs[hi]
        del cluster_members[hi]
        del cluster_scs[hi]

    return {seg.id: find(seg.id) for seg in segments}
