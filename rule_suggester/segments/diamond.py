"""Running DIAMOND all-vs-all for one product group."""

from __future__ import annotations

import subprocess
from pathlib import Path

from .params import Params


def run_diamond(group_fasta: Path, work_dir: Path, params: Params, diamond_bin: Path) -> Path:
    """makedb + blastp all-vs-all on group_fasta; returns the path to the raw outfmt-6 hits."""
    db_path = work_dir / "group"
    hits_path = work_dir / "raw_hits.tsv"

    _run(
        [
            str(diamond_bin), "makedb",
            "--in", str(group_fasta),
            "-d", str(db_path),
        ],
        step="diamond makedb",
    )
    _run(
        [
            str(diamond_bin), "blastp",
            "-q", str(group_fasta),
            "-d", str(db_path),
            "-o", str(hits_path),
            "--more-sensitive",
            "--max-target-seqs", "0",
            "--max-hsps", "0",
            "-e", str(params.evalue),
            "--outfmt", "6",
            "qseqid", "sseqid", "qstart", "qend", "sstart", "send",
            "qlen", "slen", "bitscore", "pident", "evalue",
        ],
        step="diamond blastp",
    )
    return hits_path


def _run(cmd: list[str], step: str) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{step} failed (exit {result.returncode}):\n{result.stderr}")
