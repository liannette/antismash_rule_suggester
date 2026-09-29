"""Reading phase-1 outputs (segments.tsv / groups.tsv) for phase-2 domain mapping.

Phase 2 maps onto whatever segmentation phase 1 already produced and a human reviewed, rather
than re-running DIAMOND/clustering -- so it reads the written TSVs back in, unchanged.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SegmentRecord:
    segment_id: str
    bgc_id: str
    protein_id: str
    gene_name: str
    start: int  # 0-based half-open (converted back from the 1-based-inclusive TSV)
    end: int
    length: int


@dataclass(frozen=True)
class GroupMembership:
    group_id: str
    n_bgcs: int
    bgc_id: str
    protein_id: str
    gene_name: str
    start: int  # 0-based half-open
    end: int


def load_segments_tsv(path: Path) -> list[SegmentRecord]:
    with path.open(newline="") as f:
        return [
            SegmentRecord(
                segment_id=row["segment_id"],
                bgc_id=row["bgc"],
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
                n_bgcs=int(row["n_bgcs"]),
                bgc_id=row["bgc"],
                protein_id=row["protein_id"],
                gene_name=row["gene_name"],
                start=int(row["start"]) - 1,
                end=int(row["end"]),
            )
            for row in csv.DictReader(f, delimiter="\t")
        ]
