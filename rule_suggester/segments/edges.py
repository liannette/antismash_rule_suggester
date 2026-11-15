"""Building segment-level edges from HSPs (spec step 5)."""

from __future__ import annotations

from collections import defaultdict

from .hits import LocalHit
from .io_fasta import Protein
from .params import Params
from .segmentation import Segment


def build_edges(
    proteins: list[Protein],
    segments: list[Segment],
    local_hits: dict[str, list[LocalHit]],
    selfbits: dict[str, float],
    params: Params,
) -> list[tuple[str, str, float]]:
    proteins_by_key = {p.key: p for p in proteins}
    prot_len = {p.key: p.length for p in proteins}

    segs_by_protein: dict[str, list[Segment]] = defaultdict(list)
    for seg in segments:
        segs_by_protein[seg.key].append(seg)

    def seg_selfbits(seg: Segment) -> float:
        return selfbits[seg.key] * seg.length / prot_len[seg.key]

    best_weight: dict[tuple[str, str], float] = {}

    for q_key, hits in local_hits.items():
        q_segs = segs_by_protein[q_key]
        for h in hits:
            if proteins_by_key[q_key].sc_id == proteins_by_key[h.other_key].sc_id:
                continue  # never edge segments of the same subcluster (paralog hits are already dropped)
            hsp_len = h.self_end - h.self_start
            if hsp_len <= 0:
                continue
            s_segs = segs_by_protein[h.other_key]

            for seg_q in q_segs:
                clip_start = max(seg_q.start, h.self_start)
                clip_end = min(seg_q.end, h.self_end)
                if clip_end <= clip_start:
                    continue
                recip_q = (clip_end - clip_start) / seg_q.length
                if recip_q < params.min_recip_overlap:
                    continue
                frac_of_hsp = (clip_end - clip_start) / hsp_len

                other_span = h.other_end - h.other_start
                proj_start = h.other_start + (clip_start - h.self_start) / hsp_len * other_span
                proj_end = h.other_start + (clip_end - h.self_start) / hsp_len * other_span
                if proj_end < proj_start:
                    proj_start, proj_end = proj_end, proj_start

                for seg_s in s_segs:
                    ov_start = max(proj_start, seg_s.start)
                    ov_end = min(proj_end, seg_s.end)
                    if ov_end <= ov_start:
                        continue
                    recip_s = (ov_end - ov_start) / seg_s.length
                    if recip_s < params.min_recip_overlap:
                        continue

                    denom = min(seg_selfbits(seg_q), seg_selfbits(seg_s))
                    if denom <= 0:
                        continue
                    weight = h.bitscore * frac_of_hsp / denom

                    pair = (seg_q.id, seg_s.id) if seg_q.id < seg_s.id else (seg_s.id, seg_q.id)
                    if weight > best_weight.get(pair, 0.0):
                        best_weight[pair] = weight

    return [(a, b, w) for (a, b), w in best_weight.items() if w >= params.min_w]
