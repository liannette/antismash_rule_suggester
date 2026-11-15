"""Reading phase-1 outputs (segments.tsv / groups.tsv) -- shared by any later phase that maps
onto whatever segmentation phase 1 already produced and a human reviewed, rather than re-running
DIAMOND/clustering. Both pfam/ (phase 2) and rules/ (phase 3) read these back in, unchanged.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SegmentRecord:
    segment_id: str
    sc_id: str
    protein_id: str
    gene_name: str
    start: int  # 0-based half-open (converted back from the 1-based-inclusive TSV)
    end: int
    length: int


@dataclass(frozen=True)
class GroupMembership:
    group_id: str
    n_scs: int
    sc_id: str
    protein_id: str
    gene_name: str
    start: int  # 0-based half-open
    end: int


def load_segments_tsv(path: Path) -> list[SegmentRecord]:
    with path.open(newline="") as f:
        return [
            SegmentRecord(
                segment_id=row["segment_id"],
                sc_id=row["sc"],
                protein_id=row["protein_id"],
                gene_name=row["gene_name"],
                start=int(row["start"]) - 1,
                end=int(row["end"]),
                length=int(row["length"]),
            )
            for row in csv.DictReader(f, delimiter="\t")
        ]


def load_groups_tsv(path: Path) -> list[GroupMembership]:
    with path.open(newline="") as f:
        return [
            GroupMembership(
                group_id=row["group_id"],
                n_scs=int(row["n_scs"]),
                sc_id=row["sc"],
                protein_id=row["protein_id"],
                gene_name=row["gene_name"],
                start=int(row["start"]) - 1,
                end=int(row["end"]),
            )
            for row in csv.DictReader(f, delimiter="\t")
        ]
