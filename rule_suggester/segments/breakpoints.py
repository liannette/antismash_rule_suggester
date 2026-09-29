"""Breakpoint detection and propagation (spec step 3)."""

from __future__ import annotations

from dataclasses import dataclass
from statistics import median

from .hits import LocalHit
from .io_fasta import Protein
from .params import Params


@dataclass
class BreakpointResult:
    breakpoints: dict[str, list[int]]
    iterations: int
    converged: bool


def _merge(points: set[float], length: int, tol: float) -> set[int]:
    """Collapse points within tol of each other to their median; snap points within tol of an
    end to that end. Always keeps 0 and length.
    """
    pts = sorted(points)
    clusters: list[list[float]] = []
    for p in pts:
        if clusters and p - clusters[-1][-1] <= tol:
            clusters[-1].append(p)
        else:
            clusters.append([p])

    merged = {round(median(c)) for c in clusters}
    snapped = set()
    for m in merged:
        if m <= tol:
            m = 0
        elif m >= length - tol:
            m = length
        snapped.add(m)
    snapped.add(0)
    snapped.add(length)
    return snapped


def compute_breakpoints(
    proteins: list[Protein],
    local_hits: dict[str, list[LocalHit]],
    params: Params,
) -> BreakpointResult:
    lengths = {p.key: p.length for p in proteins}

    breakpoints: dict[str, set[int]] = {}
    for p in proteins:
        raw = {0, p.length}
        for h in local_hits[p.key]:
            raw.add(h.self_start)
            raw.add(h.self_end)
        breakpoints[p.key] = _merge(raw, p.length, params.tol)

    converged = False
    iterations = 0
    for iterations in range(1, params.max_iter + 1):
        candidates: dict[str, set[float]] = {key: set() for key in breakpoints}
        for key, hits in local_hits.items():
            for h in hits:
                span = h.self_end - h.self_start
                if span <= 0:
                    continue
                for x in breakpoints[key]:
                    if h.self_start < x < h.self_end:
                        frac = (x - h.self_start) / span
                        mapped = h.other_start + frac * (h.other_end - h.other_start)
                        candidates[h.other_key].add(mapped)

        changed = False
        for key, new_pts in candidates.items():
            if not new_pts:
                continue
            merged = _merge(breakpoints[key] | new_pts, lengths[key], params.tol)
            if merged != breakpoints[key]:
                breakpoints[key] = merged
                changed = True

        if not changed:
            converged = True
            break

    return BreakpointResult(
        breakpoints={key: sorted(pts) for key, pts in breakpoints.items()},
        iterations=iterations,
        converged=converged,
    )
