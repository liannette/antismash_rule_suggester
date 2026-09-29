from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class Params:
    """All tunable parameters for the segment-grouping pipeline, in one place so a run's
    params.json is a complete record of what produced its outputs.
    """

    evalue: float = 1e-5
    min_hsp_len: int = 50
    tol: int = 25
    min_seg: int = 30
    min_recip_overlap: float = 0.5
    min_w: float = 0.2
    min_link: float = 0.2
    max_iter: int = 10

    def to_dict(self) -> dict:
        return asdict(self)
