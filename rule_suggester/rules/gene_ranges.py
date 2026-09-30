"""Estimating antiSMASH RULE CUTOFF from real genomic gene spacing, using reference GenBank
files for a product group's example BGCs (data/genbanks/<group>/<bgc_id>.gbk). Ported from the
same approach in archived/rule_suggester/gene_ranges.py.

NEIGHBOURHOOD is not estimated here -- by design it's always 0 (see build_rule.py): CONDITIONS
already requires every core group's own representative domains, so there's no extra accessory
search radius to justify widening beyond the genes CUTOFF already covers.

groups.tsv only carries protein_id and gene_name (no locus_tag), so matching against GenBank CDS
features uses those two, substring-tolerant to absorb accession-version differences between the
segment-grouping data and the GenBank qualifiers.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO

from ..phase1_io import GroupMembership


def _matches(value: str, qualifier_values: list[str]) -> bool:
    if not value:
        return False
    return any(value in q or q in value for q in qualifier_values)


def _find_cds_feature(record, protein_id: str, gene_name: str):
    """Find the GenBank CDS feature for a gene: try protein_id first, then gene name."""
    for key, value in (("protein_id", protein_id), ("gene", gene_name)):
        if not value:
            continue
        for feature in record.features:
            if feature.type == "CDS" and _matches(value, feature.qualifiers.get(key, [])):
                return feature
    return None


@dataclass(frozen=True)
class BgcSpan:
    bgc_id: str
    span_kb: float | None  # None if the .gbk file is missing or no gene could be matched
    n_genes: int
    n_matched: int
    unmatched: list[str] = field(default_factory=list)  # "gene_name (protein_id)" that couldn't be resolved


def compute_bgc_spans(memberships: list[GroupMembership], genbank_dir: Path) -> list[BgcSpan]:
    """One span per BGC in this product group: the genomic extent (kb) from the earliest gene
    start to the latest gene end, among every gene phase 1 assigned to this product group --
    core and non-core alike, since that's the real width of the known example locus, not just
    the anchor genes.
    """
    by_bgc: dict[str, list[GroupMembership]] = {}
    for m in memberships:
        by_bgc.setdefault(m.bgc_id, []).append(m)

    spans = []
    for bgc_id, genes in sorted(by_bgc.items()):
        gbk_path = Path(genbank_dir) / f"{bgc_id}.gbk"
        if not gbk_path.exists():
            spans.append(BgcSpan(bgc_id, None, len(genes), 0))
            continue

        record = SeqIO.read(gbk_path, "genbank")
        coordinates = []
        unmatched = []
        for g in genes:
            feature = _find_cds_feature(record, g.protein_id, g.gene_name)
            if feature is not None:
                coordinates.append((int(feature.location.start), int(feature.location.end)))
            else:
                unmatched.append(f"{g.gene_name or '?'} ({g.protein_id or '?'})")

        span_kb = (max(e for _, e in coordinates) - min(s for s, _ in coordinates)) / 1000 if coordinates else None
        spans.append(BgcSpan(bgc_id, span_kb, len(genes), len(coordinates), unmatched))
    return spans


def estimate_cutoff_kb(spans: list[BgcSpan], margin: float = 1.2, round_to: int = 5) -> int | None:
    """CUTOFF (kb): the largest observed genomic span across all resolvable example BGCs, scaled
    by a safety margin so the estimate isn't a knife-edge fit to exactly the widest example, then
    rounded up to the nearest round_to kb. CUTOFF needs to be at least as wide as the widest real
    example, or that BGC's own genes wouldn't be merged into one candidate cluster. None if no
    BGC's genes could be resolved at all.
    """
    resolvable = [s.span_kb for s in spans if s.span_kb is not None]
    if not resolvable:
        return None
    return math.ceil(max(resolvable) * margin / round_to) * round_to
