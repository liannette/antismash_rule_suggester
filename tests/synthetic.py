"""Helpers for building synthetic BGC FASTA directories with known ground truth.

Random sequences make homology exact and unambiguous: two sequences are homologous iff one was
derived from the other by mutate(), never by chance.
"""

from __future__ import annotations

import random
from pathlib import Path

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def random_protein(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(AMINO_ACIDS) for _ in range(length))


def mutate(sequence: str, rate: float, seed: int) -> str:
    """Substitute a `rate` fraction of positions with a different random residue."""
    rng = random.Random(seed)
    out = list(sequence)
    for i in range(len(out)):
        if rng.random() < rate:
            choices = [a for a in AMINO_ACIDS if a != out[i]]
            out[i] = rng.choice(choices)
    return "".join(out)


def write_bgc_fasta(group_dir: Path, bgc_id: str, genes: list[tuple[str, str, str]]) -> Path:
    """genes: list of (protein_id, gene_name, sequence). locus_tag mirrors protein_id."""
    path = group_dir / f"{bgc_id}.fasta"
    with path.open("w") as f:
        for protein_id, gene_name, sequence in genes:
            locus_tag = f"{bgc_id}_{protein_id}"
            f.write(f">pid:{protein_id}|gid:{gene_name}|locustag:{locus_tag}\n{sequence}\n")
    return path
