"""End-to-end orchestration: scan Pfam (and PGAP, if configured) domains once across every protein in every selected
product group, then map onto each group's own phase-1 segments/groups.

Proteins are deduplicated by (sc_id, protein_id) across groups before scanning -- a protein that
recurs in several product groups (e.g. an aminotransferase shared by Hpg and Dpg) is only searched
against the Pfam library once.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from ..phase1_io import load_groups_tsv, load_segments_tsv
from ..segments.io_fasta import Protein, discover_group_dirs, load_group_proteins
from .clans import Clan, load_pfam_clans
from .libraries import PgapLibrary, SmcogLibrary, load_pfam_names
from .mapping import build_group_of_segment, compute_group_domains, map_domains_to_segments
from .outputs import write_pfam_outputs
from .params import PfamParams
from .scan import scan_pfam_domains


@dataclass
class PfamGroupResult:
    group_name: str
    n_domain_hits: int
    n_signature_domains: int
    n_straddling: int


def run_pfam(
    fasta_dir: Path, segments_dir: Path, params: PfamParams, clans_path: Path | None = None
) -> list[PfamGroupResult]:
    clan_map: dict[str, Clan] = load_pfam_clans(clans_path) if clans_path else {}
    pgap = None
    if params.pgap_hmm_path is not None:
        if clans_path is None or params.pgap_metadata_path is None:
            raise ValueError("PGAP needs the Pfam clans file (for Pfam names) and the PGAP metadata table")
        pgap = PgapLibrary.load(params.pgap_hmm_path, params.pgap_metadata_path, load_pfam_names(clans_path))
    smcog = None
    if params.smcog_hmm_path is not None and params.smcog_profiles:
        smcog = SmcogLibrary(params.smcog_hmm_path, params.smcog_profiles)
    group_dirs = discover_group_dirs(fasta_dir)

    proteins_by_group: dict[str, list[Protein]] = {}
    unique_proteins: dict[tuple[str, str], Protein] = {}
    for group_dir in group_dirs:
        proteins = load_group_proteins(group_dir)
        proteins_by_group[group_dir.name] = proteins
        for p in proteins:
            unique_proteins.setdefault((p.sc_id, p.protein_id or p.local_key), p)

    libraries = params.hmm_path.name + (f", {pgap.hmm_path.name} ({len(pgap.models)} models)" if pgap else "")
    libraries += f", {len(smcog.profiles)} curated smCOG profile(s)" if smcog else ""
    print(f"Scanning {len(unique_proteins)} unique proteins against {libraries} ...")
    all_domains = scan_pfam_domains(list(unique_proteins.values()), params, pgap, smcog)
    domains_by_protein: dict[tuple[str, str], list] = defaultdict(list)
    for d in all_domains:
        domains_by_protein[(d.sc_id, d.protein_id)].append(d)
    n_pfam = sum(1 for d in all_domains if d.source == "Pfam")
    n_smcog = sum(1 for d in all_domains if d.source == "smCOG")
    n_pgap = len(all_domains) - n_pfam - n_smcog
    print(f"  {len(all_domains)} domain hits ({n_pfam} Pfam, {n_pgap} PGAP, {n_smcog} smCOG) "
          f"across {len(unique_proteins)} proteins")

    results = []
    for group_dir in group_dirs:
        group_name = group_dir.name
        out_dir = segments_dir / group_name
        segments_tsv, groups_tsv = out_dir / "segments.tsv", out_dir / "groups.tsv"
        if not segments_tsv.exists() or not groups_tsv.exists():
            print(f"  skipping {group_name}: no phase-1 outputs in {out_dir}", file=sys.stderr)
            continue

        segments = load_segments_tsv(segments_tsv)
        groups = load_groups_tsv(groups_tsv)
        group_of_segment = build_group_of_segment(segments, groups)
        n_scs_by_group = {g.group_id: g.n_scs for g in groups}

        group_domains_flat = [
            d
            for p in proteins_by_group[group_name]
            for d in domains_by_protein.get((p.sc_id, p.protein_id or p.local_key), [])
        ]
        assignments = map_domains_to_segments(group_domains_flat, segments, params)
        group_domain_rows = compute_group_domains(
            assignments, group_of_segment, n_scs_by_group, params, clan_map
        )

        write_pfam_outputs(
            out_dir, group_name, assignments, group_domain_rows, params,
            sorted(n_scs_by_group), clan_map,
        )

        result = PfamGroupResult(
            group_name=group_name,
            n_domain_hits=len(group_domains_flat),
            n_signature_domains=sum(1 for g in group_domain_rows if g.is_signature),
            n_straddling=sum(1 for a in assignments if a.straddling),
        )
        results.append(result)
        print(
            f"  {group_name}: {result.n_domain_hits} domain hits, "
            f"{result.n_signature_domains} signature domains, {result.n_straddling} straddling"
        )

    return results
