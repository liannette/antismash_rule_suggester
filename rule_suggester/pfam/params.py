from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class PfamParams:
    hmm_path: Path
    min_prevalence: float = 1.0  # signature-domain threshold: fraction of a group's subclusters
    straddle_threshold: float = 0.8  # a domain "belongs" to a segment if >= this fraction of it lands there
    cpus: int = 0  # 0 = let pyhmmer pick (all available cores)
    # Pfam accession (no version) -> bitscore cutoff replacing that HMM's GA; all others use GA
    cutoff_overrides: dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["hmm_path"] = str(self.hmm_path)
        return d
