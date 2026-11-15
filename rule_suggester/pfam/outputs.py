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
    all_group_ids: list[str],
    clan_map: dict[str, Clan] | None = None,
) -> None:
    clan_map = clan_map or {}
    _write_domains_tsv(out_dir / "domains.tsv", assignments, clan_map)
    _write_group_domains_tsv(out_dir / "group_domains.tsv", group_domains)
    (out_dir / "pfam_params.json").write_text(json.dumps(params.to_dict(), indent=2) + "\n")
    _write_pfam_report(out_dir / "pfam_report.md", group_name, assignments, group_domains, params, all_group_ids)


def _write_domains_tsv(path: Path, assignments: list[DomainAssignment], clan_map: dict[str, Clan]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "sc", "protein_id", "gene_name", "pfam_acc", "pfam_name", "clan_acc", "clan_name",
                "env_from", "env_to", "bitscore", "i_evalue",
                "segment_id", "frac_domain_in_segment", "frac_segment_covered", "straddling", "passes_ga",
            ]
        )
        for a in sorted(assignments, key=lambda a: (a.domain.sc_id, a.domain.protein_id, a.domain.env_from)):
            d = a.domain
            clan = clan_map.get(d.pfam_acc)
            w.writerow(
                [
                    d.sc_id, d.protein_id, d.gene_name, d.pfam_acc, d.pfam_name,
                    clan.acc if clan else "", clan.name if clan else "",
                    d.env_from + 1, d.env_to, d.bitscore, d.i_evalue,
                    a.segment_id or "", f"{a.frac_domain_in_segment:.3f}", f"{a.frac_segment_covered:.3f}",
                    a.straddling, d.passes_ga,
                ]
            )


def _write_group_domains_tsv(path: Path, group_domains: list[GroupDomain]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "group_id", "pfam_acc", "pfam_name", "clan_acc", "clan_name",
                "n_scs_with_domain", "n_scs_total", "prevalence", "is_signature", "any_straddling",
                "best_i_evalue", "is_representative",
            ]
        )
        for g in sorted(group_domains, key=lambda g: (g.group_id, -g.prevalence, g.pfam_acc)):
            w.writerow(
                [
                    g.group_id, g.pfam_acc, g.pfam_name, g.clan_acc, g.clan_name,
                    g.n_scs_with_domain, g.n_scs_total,
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
    all_group_ids: list[str],
) -> None:
    lines = [f"# Pfam domain report: {group_name}", ""]

    by_group: dict[str, list[GroupDomain]] = {}
    for g in group_domains:
        by_group.setdefault(g.group_id, []).append(g)

    no_hits: list[str] = []  # groups with zero Pfam domain hits at all -- would otherwise vanish silently
    no_signature: list[str] = []  # groups with hits, but none reaching the signature threshold

    lines.append("## Domains per group")
    lines.append("")
    lines.append(f"Signature = prevalence >= {params.min_prevalence} and not straddling.")
    lines.append(
        "Representative = the signature domain with the best (lowest) independent E-value in its "
        "clan; signature domains with no clan assignment are automatically representative."
    )
    lines.append("")
    for group_id in sorted(set(all_group_ids) | set(by_group)):
        rows = sorted(
            by_group.get(group_id, []),
            key=lambda g: (not g.is_signature, not g.is_representative, -g.prevalence, g.pfam_acc),
        )
        lines.append(f"### {group_id}")
        lines.append("")
        if rows:
            lines.append("| Pfam | name | clan | subclusters with domain | signature | representative |")
            lines.append("|---|---|---|---|---|---|")
            for g in rows:
                clan = f"{g.clan_acc} {g.clan_name}".strip() if g.clan_acc else ""
                lines.append(
                    f"| {g.pfam_acc} | {g.pfam_name} | {clan} | {g.n_scs_with_domain}/{g.n_scs_total} | "
                    f"{'yes' if g.is_signature else 'no'} | {'yes' if g.is_representative else 'no'} |"
                )
            if not any(g.is_signature for g in rows):
                no_signature.append(group_id)
                lines.append("")
                lines.append(
                    f"**WARNING: no domain in {group_id} reaches the signature threshold** -- "
                    "this group can't contribute to CONDITIONS/EXTENDERS in Step 3."
                )
        else:
            no_hits.append(group_id)
            lines.append(
                f"**WARNING: no Pfam domain hits at all in {group_id}** -- it exists in Step 1's "
                "segment grouping but has no signature/representative domain here, and can't "
                "contribute to CONDITIONS/EXTENDERS in Step 3."
            )
        lines.append("")

    lines.extend(_cutoff_override_section(assignments, params))

    straddling = [a for a in assignments if a.straddling]
    lines.append("## Straddling domains (possible misplaced breakpoint)")
    lines.append("")
    if straddling:
        lines.append("| subcluster | protein | Pfam | name | env range | best segment share |")
        lines.append("|---|---|---|---|---|---|")
        for a in sorted(straddling, key=lambda a: (a.domain.sc_id, a.domain.protein_id, a.domain.env_from)):
            d = a.domain
            lines.append(
                f"| {d.sc_id} | {d.protein_id} | {d.pfam_acc} | {d.pfam_name} | "
                f"{d.env_from + 1}-{d.env_to} | {a.frac_domain_in_segment:.2f} |"
            )
    else:
        lines.append("None.")
    lines.append("")

    lines.append("## Warnings")
    lines.append("")
    if no_hits or no_signature:
        if no_hits:
            lines.append(f"**No Pfam domain hits at all:** {', '.join(sorted(no_hits))}")
        if no_signature:
            lines.append(f"**Hits found, but none reach the signature threshold:** {', '.join(sorted(no_signature))}")
    else:
        lines.append("None.")
    lines.append("")

    path.write_text("\n".join(lines))


def _cutoff_override_section(assignments: list[DomainAssignment], params: PfamParams) -> list[str]:
    lines = ["## Cutoff overrides", ""]
    if not params.cutoff_overrides:
        return lines + ["None -- every Pfam uses its gathering threshold (GA).", ""]
    lines.append(
        "These Pfams use the bitscore cutoff below instead of GA (applied to sequence and domain score): "
        + ", ".join(f"{acc} {cut:g}" for acc, cut in sorted(params.cutoff_overrides.items()))
        + "."
    )
    lines.append("")
    below_ga = [a for a in assignments if not a.domain.passes_ga]
    if not below_ga:
        lines += ["No hit in this group depends on an override.", ""]
        return lines
    lines.append("Hits in this group detected only because of an override (below GA):")
    lines.append("")
    lines.append("| subcluster | protein | gene | Pfam | name | domain bitscore | override | segment |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for a in sorted(below_ga, key=lambda a: (a.domain.pfam_acc, a.domain.sc_id, a.domain.protein_id)):
        d = a.domain
        lines.append(
            f"| {d.sc_id} | {d.protein_id} | {d.gene_name} | {d.pfam_acc} | {d.pfam_name} | "
            f"{d.bitscore:.1f} | {params.cutoff_overrides[d.pfam_acc]:g} | {a.segment_id or '-'} |"
        )
    lines.append("")
    return lines
