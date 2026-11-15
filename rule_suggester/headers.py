"""FASTA header parsing.

export_fastas.py (subcluster_collection) writes one FASTA file per subcluster (SC####.fasta),
with headers of the form `sc:<subcluster_id>|acc:<accession>|pid:<protein_id>|gid:<gene_name>|
locustag:<locus_tag>` -- a field is blank, not omitted, when unknown; see
subcluster_collection/subcluster/fasta_export.py:cds_header. `acc` is the accession of the
sequence record the gene comes from (a MIBiG BGC accession, or in future a GenBank accession).

The same format is written into a group's concatenated multi-FASTA for DIAMOND (io_fasta.py), so
this is the only header format read anywhere -- input loading, hit parsing, segment IDs and output
FASTA headers all go through parse_header().
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ProteinHeader:
    sc_id: str
    acc: str
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


def make_header(sc_id: str, acc: str, protein_id: str, gene_name: str, locus_tag: str) -> str:
    """Build a sc/acc/pid/gid/locustag header."""
    return f"sc:{sc_id}|acc:{acc}|pid:{protein_id}|gid:{gene_name}|locustag:{locus_tag}"


def parse_header(header: str) -> ProteinHeader:
    """Parse a sc/acc/pid/gid/locustag header. Missing fields come back as ''."""
    fields = _parse_fields(header)
    return ProteinHeader(
        sc_id=fields.get("sc", ""),
        acc=fields.get("acc", ""),
        protein_id=fields.get("pid", ""),
        gene_name=fields.get("gid", ""),
        locus_tag=fields.get("locustag", ""),
    )
