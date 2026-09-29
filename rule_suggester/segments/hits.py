"""Loading and filtering DIAMOND hits (spec step 2).

Coordinates are converted from BLAST's 1-based inclusive to 0-based half-open here, once, so every
later stage works in the same coordinate system.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..headers import parse_header
from .io_fasta import Protein
from .params import Params


@dataclass(frozen=True)
class Hsp:
    """One filtered, cross-BGC HSP, coordinates 0-based half-open."""

    q_key: str
    s_key: str
    q_start: int
    q_end: int
    s_start: int
    s_end: int
    bitscore: float


@dataclass(frozen=True)
class LocalHit:
    """An HSP viewed from one of its two proteins: self_* are this protein's own coordinates,
    other_* are the partner's. Every filtered cross-BGC HSP is stored from both proteins' point of
    view, so breakpoint propagation and edge-building never have to special-case which side of the
    original DIAMOND row a protein was on.
    """

    other_key: str
    self_start: int
    self_end: int
    other_start: int
    other_end: int
    bitscore: float


@dataclass
class HitsResult:
    filtered_hits: list[Hsp]
    local_hits: dict[str, list[LocalHit]]
    selfbits: dict[str, float]
    proteins_without_self_hit: list[str]


def load_hits(hits_path: Path, proteins: list[Protein], params: Params) -> HitsResult:
    proteins_by_key = {p.key: p for p in proteins}

    self_bitscores: dict[str, list[float]] = {p.key: [] for p in proteins}
    fallback_bitscores: dict[str, list[float]] = {p.key: [] for p in proteins}
    filtered_hits: list[Hsp] = []
    local_hits: dict[str, list[LocalHit]] = {p.key: [] for p in proteins}

    with hits_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            (
                qseqid, sseqid, qstart, qend, sstart, send,
                qlen, slen, bitscore, pident, evalue,
            ) = line.split("\t")

            q = parse_header(qseqid)
            s = parse_header(sseqid)
            q_key = f"{q.bgc_id}|{q.protein_id or q.locus_tag or q.gene_name}"
            s_key = f"{s.bgc_id}|{s.protein_id or s.locus_tag or s.gene_name}"
            if q_key not in proteins_by_key or s_key not in proteins_by_key:
                raise ValueError(
                    f"hit references unknown protein ({q_key!r} or {s_key!r}) -- "
                    "hits.tsv and the FASTA it came from are out of sync"
                )

            bits = float(bitscore)
            q_start0, q_end0 = int(qstart) - 1, int(qend)
            s_start0, s_end0 = int(sstart) - 1, int(send)
            fallback_bitscores[q_key].append(bits)
            fallback_bitscores[s_key].append(bits)

            if q_key == s_key:
                self_bitscores[q_key].append(bits)
                continue

            if q.bgc_id == s.bgc_id:
                # different genes of the same BGC: paralogs never form edges.
                continue

            if (q_end0 - q_start0) < params.min_hsp_len:
                continue

            filtered_hits.append(
                Hsp(q_key, s_key, q_start0, q_end0, s_start0, s_end0, bits)
            )
            local_hits[q_key].append(
                LocalHit(s_key, q_start0, q_end0, s_start0, s_end0, bits)
            )
            local_hits[s_key].append(
                LocalHit(q_key, s_start0, s_end0, q_start0, q_end0, bits)
            )

    selfbits: dict[str, float] = {}
    proteins_without_self_hit: list[str] = []
    for key in proteins_by_key:
        if self_bitscores[key]:
            selfbits[key] = max(self_bitscores[key])
        else:
            proteins_without_self_hit.append(key)
            selfbits[key] = max(fallback_bitscores[key], default=0.0)

    return HitsResult(filtered_hits, local_hits, selfbits, proteins_without_self_hit)
