from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class PfamParams:
    hmm_path: Path
    min_prevalence: float = 1.0  # signature-domain threshold: fraction of a group's subclusters
    straddle_threshold: float = 0.8  # a domain "belongs" to a segment if >= this fraction of it lands there
    cpus: int = 0  # 0 = let pyhmmer pick (all available cores)
    # model accession (no version) -> bitscore cutoff replacing that HMM's GA; all others use GA
    cutoff_overrides: dict[str, float] = field(default_factory=dict)
    pgap_hmm_path: Path | None = None  # NCBI PGAP library (hmm_PGAP.HMM.tgz); None = Pfam only
    pgap_metadata_path: Path | None = None  # hmm_PGAP.tsv
    # a PGAP representative replaces a Pfam signature domain if it covers >= this share of the
    # Pfam domain's envelope (in a majority of the group's subclusters); also used between PGAP models
    replace_overlap: float = 0.5
    # PGAP libraries whose signature models can become representatives (and so replace Pfam);
    # PRK models are full-length, mostly primary-metabolism clusters with strict cutoffs, so they're
    # reported but not used in rules
    representative_pgap_sources: tuple[str, ...] = ("TIGRFAM", "NCBIFAM")

    def to_dict(self) -> dict:
        d = asdict(self)
        for key in ("hmm_path", "pgap_hmm_path", "pgap_metadata_path"):
            d[key] = str(d[key]) if d[key] is not None else None
        return d
