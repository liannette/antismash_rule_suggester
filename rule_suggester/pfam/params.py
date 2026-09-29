from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class PfamParams:
    hmm_path: Path
    min_prevalence: float = 0.8  # signature-domain threshold: fraction of a group's BGCs
    straddle_threshold: float = 0.8  # a domain "belongs" to a segment if >= this fraction of it lands there
    cpus: int = 0  # 0 = let pyhmmer pick (all available cores)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["hmm_path"] = str(self.hmm_path)
        return d
