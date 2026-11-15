"""Estimating antiSMASH RULE CUTOFF from real genomic gene spacing, using reference GenBank
files for a product group's example subclusters (data/genbanks/<group>/<sc_id>.gbk). The CUTOFF is 
estimated as an anchor radius: antiSMASH looks CUTOFF kb left and right of each anchor CDS.

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
class ScSpan:
    sc_id: str
    span_kb: float | None  # None if the .gbk file is missing or no gene could be matched
    n_genes: int
    n_matched: int
    unmatched: list[str] = field(default_factory=list)  # "gene_name (protein_id)" that couldn't be resolved
    radius_kb: float | None = None  # None if no CONDITIONS gene could be matched


def _edge_distance(a: tuple[int, int], b: tuple[int, int]) -> int:
    """Gap (bp) between two genes' nearest edges; 0 if they touch or overlap."""
    return max(0, b[0] - a[1], a[0] - b[1])


def _anchor_radius(coordinates: list[tuple[int, int]]) -> int:
    """Smallest reach (bp) any one gene needs, measured outward from its own edges, to cover every
    other gene -- i.e. the best anchor's worst-case distance. antiSMASH evaluates CONDITIONS on the
    genes within CUTOFF of an anchor CDS on either side, so this is the tightest CUTOFF under which
    some anchor sees all of them at once.
    """
    return min(max(_edge_distance(a, g) for g in coordinates) for a in coordinates)


def compute_sc_spans(
    memberships: list[GroupMembership], genbank_dir: Path, condition_groups: set[str]
) -> list[ScSpan]:
    """One entry per subcluster in this product group, with two extents:

    - span_kb: earliest gene start to latest gene end, among every gene phase 1 assigned to this
      product group (core and non-core alike) -- the real width of the known example locus, kept
      for reference.
    - radius_kb: the anchor radius (see _anchor_radius) over only the genes in condition_groups,
      i.e. the groups CONDITIONS requires. That's what CUTOFF is estimated from: genes outside
      CONDITIONS don't need to be within CUTOFF for the rule to fire.
    """
    by_sc: dict[str, list[GroupMembership]] = {}
    for m in memberships:
        by_sc.setdefault(m.sc_id, []).append(m)

    spans = []
    for sc_id, genes in sorted(by_sc.items()):
        gbk_path = Path(genbank_dir) / f"{sc_id}.gbk"
        if not gbk_path.exists():
            spans.append(ScSpan(sc_id, None, len(genes), 0))
            continue

        record = SeqIO.read(gbk_path, "genbank")
        coordinates = []
        condition_coordinates = []
        unmatched = []
        for g in genes:
            feature = _find_cds_feature(record, g.protein_id, g.gene_name)
            if feature is not None:
                location = (int(feature.location.start), int(feature.location.end))
                coordinates.append(location)
                if g.group_id in condition_groups:
                    condition_coordinates.append(location)
            else:
                unmatched.append(f"{g.gene_name or '?'} ({g.protein_id or '?'})")

        span_kb = (max(e for _, e in coordinates) - min(s for s, _ in coordinates)) / 1000 if coordinates else None
        radius_kb = _anchor_radius(condition_coordinates) / 1000 if condition_coordinates else None
        spans.append(ScSpan(sc_id, span_kb, len(genes), len(coordinates), unmatched, radius_kb))
    return spans


def estimate_cutoff_kb(spans: list[ScSpan], round_to: int = 5) -> int | None:
    """CUTOFF (kb): the largest anchor radius across all resolvable example subclusters, rounded up to the
    nearest round_to kb. CUTOFF needs to be at least every example's radius, or that subcluster's
    CONDITIONS genes would never all fall within CUTOFF of a single anchor. At least round_to kb:
    genes that touch in every example give radius 0, but CUTOFF 0 would miss the same genes in a
    genome with even a small gap between them. None if no subcluster's CONDITIONS genes could be
    resolved at all.
    """
    resolvable = [s.radius_kb for s in spans if s.radius_kb is not None]
    if not resolvable:
        return None
    return max(round_to, math.ceil(max(resolvable) / round_to) * round_to)
