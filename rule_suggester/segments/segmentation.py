"""Turning breakpoints into segments (spec step 4)."""

from __future__ import annotations

from dataclasses import dataclass

from .io_fasta import Protein
from .params import Params


@dataclass(frozen=True)
class Segment:
    id: str
    key: str  # parent protein key
    bgc_id: str
    protein_id: str
    gene_name: str
    locus_tag: str
    start: int  # 0-based half-open
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


def _merge_short(breakpoints: list[int], min_seg: int) -> list[int]:
    """Repeatedly merge the shortest segment into whichever of its neighbours is shorter, until
    no segment is shorter than min_seg or only one segment remains.
    """
    pts = list(breakpoints)
    while len(pts) > 2:
        lengths = [pts[i + 1] - pts[i] for i in range(len(pts) - 1)]
        shortest_len = min(lengths)
        if shortest_len >= min_seg:
            break
        i = lengths.index(shortest_len)
        left_len = lengths[i - 1] if i - 1 >= 0 else None
        right_len = lengths[i + 1] if i + 1 < len(lengths) else None
        if right_len is None or (left_len is not None and left_len <= right_len):
            del pts[i]  # merge into the left neighbour
        else:
            del pts[i + 1]  # merge into the right neighbour
    return pts


def build_segments(
    proteins: list[Protein],
    breakpoints: dict[str, list[int]],
    params: Params,
) -> list[Segment]:
    segments: list[Segment] = []
    for p in proteins:
        pts = _merge_short(breakpoints[p.key], params.min_seg)
        pid = p.protein_id or p.local_key
        for start, end in zip(pts, pts[1:]):
            seg_id = f"{p.bgc_id}|{pid}|{start + 1}-{end}"
            segments.append(
                Segment(
                    id=seg_id,
                    key=p.key,
                    bgc_id=p.bgc_id,
                    protein_id=p.protein_id,
                    gene_name=p.gene_name,
                    locus_tag=p.locus_tag,
                    start=start,
                    end=end,
                )
            )
    return segments
