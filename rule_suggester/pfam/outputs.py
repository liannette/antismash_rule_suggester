"""Writing Pfam domain outputs for one product group (spec Phase 2)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .mapping import DomainAssignment, GroupDomain
from .params import PfamParams


def write_pfam_outputs(
    out_dir: Path,
    group_name: str,
    assignments: list[DomainAssignment],
    group_domains: list[GroupDomain],
    params: PfamParams,
) -> None:
    _write_domains_tsv(out_dir / "domains.tsv", assignments)
    _write_group_domains_tsv(out_dir / "group_domains.tsv", group_domains)
    (out_dir / "pfam_params.json").write_text(json.dumps(params.to_dict(), indent=2) + "\n")
    _write_pfam_report(out_dir / "pfam_report.md", group_name, assignments, group_domains, params)


def _write_domains_tsv(path: Path, assignments: list[DomainAssignment]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "bgc", "protein_id", "gene_name", "pfam_acc", "pfam_name",
                "env_from", "env_to", "bitscore", "i_evalue",
                "segment_id", "frac_domain_in_segment", "frac_segment_covered", "straddling",
            ]
        )
        for a in sorted(assignments, key=lambda a: (a.domain.bgc_id, a.domain.protein_id, a.domain.env_from)):
            d = a.domain
            w.writerow(
                [
                    d.bgc_id, d.protein_id, d.gene_name, d.pfam_acc, d.pfam_name,
                    d.env_from + 1, d.env_to, d.bitscore, d.i_evalue,
                    a.segment_id or "", f"{a.frac_domain_in_segment:.3f}", f"{a.frac_segment_covered:.3f}",
                    a.straddling,
                ]
            )


def _write_group_domains_tsv(path: Path, group_domains: list[GroupDomain]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "group_id", "pfam_acc", "pfam_name", "n_bgcs_with_domain", "n_bgcs_total",
                "prevalence", "is_signature", "any_straddling",
            ]
        )
        for g in sorted(group_domains, key=lambda g: (g.group_id, -g.prevalence, g.pfam_acc)):
            w.writerow(
                [
                    g.group_id, g.pfam_acc, g.pfam_name, g.n_bgcs_with_domain, g.n_bgcs_total,
                    f"{g.prevalence:.3f}", g.is_signature, g.any_straddling,
                ]
            )


def _write_pfam_report(
    path: Path,
    group_name: str,
    assignments: list[DomainAssignment],
    group_domains: list[GroupDomain],
    params: PfamParams,
) -> None:
    lines = [f"# Pfam domain report: {group_name}", ""]

    by_group: dict[str, list[GroupDomain]] = {}
    for g in group_domains:
        by_group.setdefault(g.group_id, []).append(g)

    lines.append("## Signature domains per group")
    lines.append("")
    lines.append(f"(prevalence >= {params.min_prevalence}, not straddling)")
    lines.append("")
    for group_id in sorted(by_group):
        sig = sorted((g for g in by_group[group_id] if g.is_signature), key=lambda g: -g.prevalence)
        lines.append(f"### {group_id}")
        lines.append("")
        if sig:
            lines.append("| Pfam | name | BGCs with domain |")
            lines.append("|---|---|---|")
            for g in sig:
                lines.append(f"| {g.pfam_acc} | {g.pfam_name} | {g.n_bgcs_with_domain}/{g.n_bgcs_total} |")
        else:
            lines.append("None.")
        lines.append("")

    straddling = [a for a in assignments if a.straddling]
    lines.append("## Straddling domains (possible misplaced breakpoint)")
    lines.append("")
    if straddling:
        lines.append("| BGC | protein | Pfam | name | env range | best segment share |")
        lines.append("|---|---|---|---|---|---|")
        for a in sorted(straddling, key=lambda a: (a.domain.bgc_id, a.domain.protein_id, a.domain.env_from)):
            d = a.domain
            lines.append(
                f"| {d.bgc_id} | {d.protein_id} | {d.pfam_acc} | {d.pfam_name} | "
                f"{d.env_from + 1}-{d.env_to} | {a.frac_domain_in_segment:.2f} |"
            )
    else:
        lines.append("None.")
    lines.append("")

    path.write_text("\n".join(lines))
