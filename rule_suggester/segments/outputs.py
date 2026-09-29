"""Writing every output file for one product group (spec 'Outputs' section)."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from .absorb import SegmentMerge
from .breakpoints import BreakpointResult
from .grouping import Group
from .hits import HitsResult
from .io_fasta import Protein
from .params import Params
from .segmentation import Segment


def write_outputs(
    out_dir: Path,
    group_name: str,
    proteins: list[Protein],
    hits: HitsResult,
    breakpoint_result: BreakpointResult,
    segments: list[Segment],
    edges: list[tuple[str, str, float]],
    groups: list[Group],
    invariant_violations: list[str],
    absorptions: list[SegmentMerge],
    params: Params,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    _write_hits_tsv(out_dir / "hits.tsv", hits.filtered_hits)
    _write_segments_tsv(out_dir / "segments.tsv", segments)
    _write_edges_tsv(out_dir / "edges.tsv", edges)
    _write_groups_tsv(out_dir / "groups.tsv", groups)
    _write_group_fastas(out_dir / "groups", groups, proteins)
    _write_presence_tsv(out_dir / "presence.tsv", proteins, groups)
    (out_dir / "params.json").write_text(json.dumps(params.to_dict(), indent=2) + "\n")
    _write_report(
        out_dir / "report.md",
        group_name,
        proteins,
        hits,
        breakpoint_result,
        segments,
        groups,
        invariant_violations,
        absorptions,
    )


def _write_hits_tsv(path: Path, hits) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["q_key", "s_key", "q_start", "q_end", "s_start", "s_end", "bitscore"])
        for h in hits:
            w.writerow([h.q_key, h.s_key, h.q_start + 1, h.q_end, h.s_start + 1, h.s_end, h.bitscore])


def _segments_by_protein(segments: list[Segment]) -> dict[str, list[Segment]]:
    by_protein: dict[str, list[Segment]] = defaultdict(list)
    for seg in segments:
        by_protein[seg.key].append(seg)
    return by_protein


def _write_segments_tsv(path: Path, segments: list[Segment]) -> None:
    by_protein = _segments_by_protein(segments)
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            ["segment_id", "bgc", "protein_id", "gene_name", "start", "end", "length", "n_segments_in_protein"]
        )
        for seg in sorted(segments, key=lambda s: (s.bgc_id, s.key, s.start)):
            n_in_protein = len(by_protein[seg.key])
            w.writerow(
                [seg.id, seg.bgc_id, seg.protein_id, seg.gene_name, seg.start + 1, seg.end, seg.length, n_in_protein]
            )


def _write_edges_tsv(path: Path, edges: list[tuple[str, str, float]]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["seg_a", "seg_b", "weight"])
        for a, b, weight in sorted(edges, key=lambda e: (-e[2], e[0], e[1])):
            w.writerow([a, b, weight])


def _write_groups_tsv(path: Path, groups: list[Group]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["group_id", "n_bgcs", "bgc", "protein_id", "gene_name", "start", "end", "seg_len"])
        for group in groups:
            for seg in group.segments:
                w.writerow(
                    [group.id, group.n_bgcs, seg.bgc_id, seg.protein_id, seg.gene_name, seg.start + 1, seg.end, seg.length]
                )


def _write_group_fastas(groups_dir: Path, groups: list[Group], proteins: list[Protein]) -> None:
    groups_dir.mkdir(parents=True, exist_ok=True)
    proteins_by_key = {p.key: p for p in proteins}
    for group in groups:
        with (groups_dir / f"{group.id}.faa").open("w") as f:
            for seg in group.segments:
                protein = proteins_by_key[seg.key]
                seq = protein.sequence[seg.start:seg.end]
                header = f"{seg.bgc_id}|{seg.protein_id or protein.local_key}|{seg.gene_name}|{seg.start + 1}-{seg.end}"
                f.write(f">{header}\n{seq}\n")


def _write_presence_tsv(path: Path, proteins: list[Protein], groups: list[Group]) -> None:
    bgcs = sorted({p.bgc_id for p in proteins})
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["bgc"] + [g.id for g in groups])
        for bgc in bgcs:
            row = [bgc]
            for group in groups:
                present = any(seg.bgc_id == bgc for seg in group.segments)
                row.append(1 if present else 0)
            w.writerow(row)


def _write_report(
    path: Path,
    group_name: str,
    proteins: list[Protein],
    hits: HitsResult,
    breakpoint_result: BreakpointResult,
    segments: list[Segment],
    groups: list[Group],
    invariant_violations: list[str],
    absorptions: list[SegmentMerge],
) -> None:
    bgcs = sorted({p.bgc_id for p in proteins})
    total_bgcs = len(bgcs)
    by_protein = _segments_by_protein(segments)
    core_groups = [g for g in groups if g.n_bgcs == total_bgcs]
    singleton_groups = [g for g in groups if len(g.segments) == 1]
    split_proteins = [p for p in proteins if len(by_protein[p.key]) > 1]
    group_of_segment = {seg.id: g.id for g in groups for seg in g.segments}

    lines: list[str] = []
    lines.append(f"# Segment grouping report: {group_name}")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- BGCs: {total_bgcs}")
    lines.append(f"- Proteins: {len(proteins)}")
    lines.append(f"- Segments: {len(segments)}")
    lines.append(f"- Groups: {len(groups)}")
    lines.append(f"- Core groups (present in all BGCs): {len(core_groups)}")
    lines.append(f"- Singleton groups: {len(singleton_groups)}")
    lines.append(f"- Breakpoint propagation: {breakpoint_result.iterations} iteration(s), "
                 f"{'converged' if breakpoint_result.converged else 'did NOT converge'}")
    lines.append("")

    lines.append("## Invariant checks")
    lines.append("")
    if invariant_violations:
        lines.append(f"**{len(invariant_violations)} violation(s):**")
        lines.append("")
        for v in invariant_violations:
            lines.append(f"- {v}")
    else:
        lines.append("All invariants passed.")
    lines.append("")

    lines.append("## Groups")
    lines.append("")
    for group in groups:
        core = " (core)" if group.n_bgcs == total_bgcs else ""
        lines.append(f"### {group.id} -- {group.n_bgcs}/{total_bgcs} BGCs{core}")
        lines.append("")
        lines.append("| BGC | gene name | protein ID | range | length |")
        lines.append("|---|---|---|---|---|")
        for seg in sorted(group.segments, key=lambda s: s.bgc_id):
            lines.append(
                f"| {seg.bgc_id} | {seg.gene_name} | {seg.protein_id} | "
                f"{seg.start + 1}-{seg.end} | {seg.length} |"
            )
        lines.append("")

    lines.append("## Uninformative breakpoints merged")
    lines.append("")
    lines.append(
        "Two adjacent same-protein segments were merged (boundary extended, breakpoint dropped) "
        "whenever their adjacency covered every member of at least one of their two groups -- i.e. "
        "no BGC had one piece without the other, so the split didn't distinguish anything."
    )
    lines.append("")
    if absorptions:
        lines.append("| protein | merged segments | into segment |")
        lines.append("|---|---|---|")
        proteins_by_key = {p.key: p for p in proteins}
        for a in absorptions:
            p = proteins_by_key[a.protein_key]
            lines.append(f"| {p.protein_id or p.local_key} ({p.bgc_id}) | {a.from_ids} | {a.into_id} |")
    else:
        lines.append("None.")
    lines.append("")

    lines.append("## Proteins split into multiple segments")
    lines.append("")
    if split_proteins:
        lines.append("| protein | BGC | n segments | segments (range -> group) |")
        lines.append("|---|---|---|---|")
        for p in split_proteins:
            segs = sorted(by_protein[p.key], key=lambda s: s.start)
            detail = "; ".join(
                f"{s.start + 1}-{s.end} -> {group_of_segment.get(s.id, '?')}" for s in segs
            )
            lines.append(f"| {p.protein_id or p.local_key} | {p.bgc_id} | {len(segs)} | {detail} |")
    else:
        lines.append("None.")
    lines.append("")

    lines.append("## Flags")
    lines.append("")
    flagged = False
    noisy = [p for p in proteins if len(by_protein[p.key]) > 4]
    if noisy:
        flagged = True
        lines.append("**Proteins with >4 segments (likely breakpoint noise):**")
        for p in noisy:
            lines.append(f"- {p.protein_id or p.local_key} ({p.bgc_id}): {len(by_protein[p.key])} segments")
    if not breakpoint_result.converged:
        flagged = True
        lines.append(
            f"**Breakpoint propagation did not converge** after {breakpoint_result.iterations} iterations."
        )
    if hits.proteins_without_self_hit:
        flagged = True
        lines.append("**Proteins with no self hit (fallback bitscore used):**")
        for key in hits.proteins_without_self_hit:
            lines.append(f"- {key}")
    if not flagged:
        lines.append("None.")
    lines.append("")

    path.write_text("\n".join(lines))
