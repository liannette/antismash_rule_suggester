"""End-to-end orchestration for one product group, and for a --fasta-dir of many."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path

from .absorb import merge_uninformative_breakpoints
from .breakpoints import compute_breakpoints
from .clustering import cluster_segments
from .diamond import run_diamond
from .edges import build_edges
from .grouping import assign_groups
from .hits import load_hits
from .invariants import check_invariants
from .io_fasta import discover_group_dirs, load_group_proteins, write_group_fasta
from .outputs import write_outputs
from .params import Params
from .segmentation import build_segments


@dataclass
class GroupRunResult:
    group_name: str
    out_dir: Path
    n_proteins: int
    n_segments: int
    n_groups: int
    invariant_violations: list[str]


def run_group(group_dir: Path, out_dir: Path, params: Params, diamond_bin: Path) -> GroupRunResult:
    group_name = group_dir.name
    proteins = load_group_proteins(group_dir)
    total_bgcs = len({p.bgc_id for p in proteins})

    with tempfile.TemporaryDirectory(prefix=f"diamond_{group_name}_") as tmp:
        tmp_path = Path(tmp)
        group_fasta = write_group_fasta(proteins, tmp_path / "group.faa")
        raw_hits_path = run_diamond(group_fasta, tmp_path, params, diamond_bin)
        hits = load_hits(raw_hits_path, proteins, params)

    breakpoint_result = compute_breakpoints(proteins, hits.local_hits, params)
    segments = build_segments(proteins, breakpoint_result.breakpoints, params)
    edges = build_edges(proteins, segments, hits.local_hits, hits.selfbits, params)
    cluster_of = cluster_segments(segments, edges, params)
    segments, cluster_of, _ = merge_uninformative_breakpoints(proteins, segments, cluster_of)
    groups = assign_groups(segments, cluster_of, proteins, total_bgcs)
    violations = check_invariants(proteins, segments, groups)

    write_outputs(
        out_dir, group_name, proteins, hits, breakpoint_result, segments, edges, groups, violations, params
    )

    return GroupRunResult(
        group_name=group_name,
        out_dir=out_dir,
        n_proteins=len(proteins),
        n_segments=len(segments),
        n_groups=len(groups),
        invariant_violations=violations,
    )


def run_all(fasta_dir: Path, out_dir: Path, params: Params, diamond_bin: Path) -> list[GroupRunResult]:
    results = []
    for group_dir in discover_group_dirs(fasta_dir):
        result = run_group(group_dir, out_dir / group_dir.name, params, diamond_bin)
        results.append(result)
    return results
