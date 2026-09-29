"""FASTA header parsing.

export_fastas.py (subcluster_collection) writes one FASTA file per BGC, with headers of the form
`pid:<protein_id>|gid:<gene_name>|locustag:<locus_tag>` (a field is blank, not omitted, when
unknown) -- see subcluster_collection/subcluster/fasta_export.py:cds_header. That header alone
doesn't carry a BGC id; it's implied by which file it came from.

When a group's per-BGC files are concatenated into one multi-FASTA for DIAMOND (io_fasta.py), the
BGC id is folded into a new header `bgc:<bgc_id>|pid:...|gid:...|locustag:...`. From that point on
this is the only header format read anywhere -- hit parsing, segment IDs, and output FASTA headers
all call parse_header() on it rather than re-deriving the BGC id from a filename.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ProteinHeader:
    bgc_id: str
    protein_id: str
    gene_name: str
    locus_tag: str


def _parse_fields(header: str) -> dict[str, str]:
    header = header.lstrip(">").strip()
    fields: dict[str, str] = {}
    for part in header.split("|"):
        key, _, value = part.partition(":")
        fields[key] = value
    return fields


def parse_source_header(header: str) -> tuple[str, str, str]:
    """Parse a pid/gid/locustag header as written by export_fastas.py (no BGC id)."""
    fields = _parse_fields(header)
    return fields.get("pid", ""), fields.get("gid", ""), fields.get("locustag", "")


def make_combined_header(bgc_id: str, protein_id: str, gene_name: str, locus_tag: str) -> str:
    """Build the bgc/pid/gid/locustag header written into a group's concatenated multi-FASTA."""
    return f"bgc:{bgc_id}|pid:{protein_id}|gid:{gene_name}|locustag:{locus_tag}"


def parse_header(header: str) -> ProteinHeader:
    """Parse a combined bgc/pid/gid/locustag header -- the single parser all later code uses."""
    fields = _parse_fields(header)
    return ProteinHeader(
        bgc_id=fields.get("bgc", ""),
        protein_id=fields.get("pid", ""),
        gene_name=fields.get("gid", ""),
        locus_tag=fields.get("locustag", ""),
    )
