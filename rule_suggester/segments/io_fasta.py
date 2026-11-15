"""Resolving --fasta-dir and loading a product group's proteins from it.

data/fasta/<group>/<SC>.fasta: one directory per product group (building block), each holding one
protein FASTA per member subcluster (see docs/segment_grouping_spec.md). --fasta-dir may point at a single such group
directory, or at a parent directory holding several of them -- detected from contents, never from
a flag.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from ..headers import make_header, parse_header

_FASTA_GLOBS = ("*.fasta", "*.faa", "*.fa")


@dataclass(frozen=True)
class Protein:
    key: str  # unique within the group: "<sc_id>|<local_key>"
    sc_id: str
    acc: str  # accession of the sequence record the gene comes from (MIBiG/GenBank)
    protein_id: str
    gene_name: str
    locus_tag: str
    local_key: str  # protein_id, or locus_tag/gene_name fallback when protein_id is blank
    order: int  # 0-based position within its subcluster's FASTA file, for group-ordering heuristics
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def header(self) -> str:
        return make_header(self.sc_id, self.acc, self.protein_id, self.gene_name, self.locus_tag)


def _group_fasta_files(directory: Path) -> list[Path]:
    return sorted({p for pattern in _FASTA_GLOBS for p in directory.glob(pattern)})


def is_group_dir(path: Path) -> bool:
    """True if path directly holds FASTA files (a single product-group directory)."""
    return len(_group_fasta_files(path)) > 0


def discover_group_dirs(fasta_dir: Path) -> list[Path]:
    """Resolve --fasta-dir to the group directories to process, in a stable (sorted) order."""
    fasta_dir = Path(fasta_dir)
    if not fasta_dir.is_dir():
        raise ValueError(f"--fasta-dir {fasta_dir} is not a directory")
    if is_group_dir(fasta_dir):
        return [fasta_dir]
    group_dirs = sorted(p for p in fasta_dir.iterdir() if p.is_dir() and is_group_dir(p))
    if not group_dirs:
        raise ValueError(
            f"{fasta_dir} contains neither FASTA files nor group subdirectories containing them"
        )
    return group_dirs


def load_group_proteins(group_dir: Path) -> list[Protein]:
    """Load every protein in a group directory's per-subcluster FASTA files, in file/record order.
    The subcluster id comes from the header's sc: field, falling back to the file name (SC####)
    for headers without one."""
    proteins: list[Protein] = []
    seen_keys: set[str] = set()
    for fasta_path in _group_fasta_files(group_dir):
        for order, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
            h = parse_header(record.description)
            sc_id = h.sc_id or fasta_path.stem
            protein_id, gene_name, locus_tag = h.protein_id, h.gene_name, h.locus_tag
            local_key = protein_id or locus_tag or gene_name
            if not local_key:
                raise ValueError(
                    f"{fasta_path}: header {record.description!r} has no protein_id, "
                    "locus_tag, or gene_name to identify it by"
                )
            key = f"{sc_id}|{local_key}"
            if key in seen_keys:
                raise ValueError(
                    f"{fasta_path}: duplicate protein identifier {local_key!r} within group "
                    f"{group_dir.name}"
                )
            seen_keys.add(key)
            proteins.append(
                Protein(
                    key=key,
                    sc_id=sc_id,
                    acc=h.acc,
                    protein_id=protein_id,
                    gene_name=gene_name,
                    locus_tag=locus_tag,
                    local_key=local_key,
                    order=order,
                    sequence=str(record.seq),
                )
            )
    if not proteins:
        raise ValueError(f"{group_dir}: no protein sequences found")
    return proteins


def write_group_fasta(proteins: list[Protein], out_path: Path) -> Path:
    """Write the concatenated, combined-header multi-FASTA that DIAMOND runs against."""
    with out_path.open("w") as f:
        for p in proteins:
            f.write(f">{p.header}\n{p.sequence}\n")
    return out_path
