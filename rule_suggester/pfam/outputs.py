"""Writing Pfam domain outputs for one product group (spec Phase 2)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .clans import Clan
from .mapping import DomainAssignment, GroupDomain
from .params import PfamParams


def write_pfam_outputs(
    out_dir: Path,
    group_name: str,
    assignments: list[DomainAssignment],
    group_domains: list[GroupDomain],
    params: PfamParams,
    clan_map: dict[str, Clan] | None = None,
) -> None:
    clan_map = clan_map or {}
    _write_domains_tsv(out_dir / "domains.tsv", assignments, clan_map)
    _write_group_domains_tsv(out_dir / "group_domains.tsv", group_domains)
    (out_dir / "pfam_params.json").write_text(json.dumps(params.to_dict(), indent=2) + "\n")
    _write_pfam_report(out_dir / "pfam_report.md", group_name, assignments, group_domains, params)


def _write_domains_tsv(path: Path, assignments: list[DomainAssignment], clan_map: dict[str, Clan]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "bgc", "protein_id", "gene_name", "pfam_acc", "pfam_name", "clan_acc", "clan_name",
                "env_from", "env_to", "bitscore", "i_evalue",
                "segment_id", "frac_domain_in_segment", "frac_segment_covered", "straddling",
            ]
        )
        for a in sorted(assignments, key=lambda a: (a.domain.bgc_id, a.domain.protein_id, a.domain.env_from)):
            d = a.domain
            clan = clan_map.get(d.pfam_acc)
            w.writerow(
                [
                    d.bgc_id, d.protein_id, d.gene_name, d.pfam_acc, d.pfam_name,
                    clan.acc if clan else "", clan.name if clan else "",
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
                "group_id", "pfam_acc", "pfam_name", "clan_acc", "clan_name",
                "n_bgcs_with_domain", "n_bgcs_total", "prevalence", "is_signature", "any_straddling",
                "best_i_evalue", "is_representative",
            ]
        )
        for g in sorted(group_domains, key=lambda g: (g.group_id, -g.prevalence, g.pfam_acc)):
            w.writerow(
                [
                    g.group_id, g.pfam_acc, g.pfam_name, g.clan_acc, g.clan_name,
                    g.n_bgcs_with_domain, g.n_bgcs_total,
                    f"{g.prevalence:.3f}", g.is_signature, g.any_straddling,
                    g.best_i_evalue, g.is_representative,
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

    lines.append("## Domains per group")
    lines.append("")
    lines.append(f"Signature = prevalence >= {params.min_prevalence} and not straddling.")
    lines.append(
        "Representative = the signature domain with the best (lowest) independent E-value in its "
        "clan; signature domains with no clan assignment are automatically representative."
    )
    lines.append("")
    for group_id in sorted(by_group):
        rows = sorted(
            by_group[group_id],
            key=lambda g: (not g.is_signature, not g.is_representative, -g.prevalence, g.pfam_acc),
        )
        lines.append(f"### {group_id}")
        lines.append("")
        if rows:
            lines.append("| Pfam | name | clan | BGCs with domain | signature | representative |")
            lines.append("|---|---|---|---|---|---|")
            for g in rows:
                clan = f"{g.clan_acc} {g.clan_name}".strip() if g.clan_acc else ""
                lines.append(
                    f"| {g.pfam_acc} | {g.pfam_name} | {clan} | {g.n_bgcs_with_domain}/{g.n_bgcs_total} | "
                    f"{'yes' if g.is_signature else 'no'} | {'yes' if g.is_representative else 'no'} |"
                )
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
