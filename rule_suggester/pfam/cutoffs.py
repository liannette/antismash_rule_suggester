"""Per-Pfam bitscore cutoffs that replace Pfam's gathering threshold (GA) for domain detection.

Overrides file: tab-separated, header `pfam_acc  cutoff  note`, `#` lines are comments.
`pfam_acc` is the accession without version (PF00550). `cutoff` is a bitscore applied to both the
sequence and the domain score, i.e. it replaces both GA values of that HMM. Every other Pfam keeps
its GA.
"""

from __future__ import annotations

import csv
from pathlib import Path


def load_cutoff_overrides(path: Path) -> dict[str, float]:
    overrides: dict[str, float] = {}
    with Path(path).open(newline="") as f:
        rows = csv.DictReader((line for line in f if line.strip() and not line.startswith("#")), delimiter="\t")
        if rows.fieldnames is None or not {"pfam_acc", "cutoff"} <= set(rows.fieldnames):
            raise ValueError(f"{path}: needs a header with at least 'pfam_acc' and 'cutoff' columns")
        for row in rows:
            acc = row["pfam_acc"].strip().split(".")[0]
            if acc in overrides:
                raise ValueError(f"{path}: {acc} listed more than once")
            try:
                overrides[acc] = float(row["cutoff"])
            except ValueError:
                raise ValueError(f"{path}: cutoff for {acc} is not a number: {row['cutoff']!r}") from None
    return overrides
