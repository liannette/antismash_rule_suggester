#!/usr/bin/env python3
"""List every gene (segment) with no "normal" HMM hit after step 2 (scripts/pfam_domains.py): no
model from any searched library (Pfam, TIGRFAM, NCBIFAM, PRK) at its gathering threshold (GA).

Hits that only pass because of data/pfam/cutoff_overrides.tsv (domains.tsv passes_ga = False) don't
count as normal. For each such gene, Pfam and PGAP are searched again without cutoffs on its
full-length protein, and the hits overlapping its segment are reported with their GA, plus whether
a cutoff override now recovers the gene.

One row per gene, written to output/no_ga_hmm_hits.tsv by default. Columns:
  group_has_ga_signature  whether the group has a signature model at GA anyway (from other genes)
  group_ga_models         models hitting the group's genes at GA, with subcluster counts
  status                  recovered_by_override / weak_hit_only / no_hit
  override_hit            the best hit passing its override cutoff (name [library] score/cutoff)
  best_weak_*             the hit overlapping the segment that comes closest to its model's GA:
                          ranked by domain bitscore / domain GA (best_weak_frac_ga), since GAs
                          differ a lot between libraries (~20-30 Pfam, 100s PGAP)
  other_weak_hits         up to three further hits in the same order, "name [library] score/GA"
"""

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pyhmmer
from Bio import SeqIO
from pyhmmer.easel import Alphabet, TextSequence

from rule_suggester.pfam.cutoffs import load_cutoff_overrides
from rule_suggester.pfam.libraries import PgapLibrary, iter_hmm_file, load_pfam_names
from rule_suggester.segments.io_fasta import load_group_proteins

ROOT = Path(__file__).resolve().parent.parent
MIN_OVERLAP = 0.5  # share of a hit's envelope that must fall in the segment for it to count there
MAX_I_EVALUE = 1.0  # relaxed search: report domains up to this i-Evalue (Z = number of models in the library)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--segments-dir", type=Path, default=ROOT / "output" / "segments")
    parser.add_argument("--fasta-root", type=Path, default=ROOT / "data" / "fasta")
    parser.add_argument("--genbank-root", type=Path, default=ROOT / "data" / "genbanks")
    parser.add_argument("--pfam-hmm", type=Path, default=ROOT / "data" / "pfam" / "Pfam-A.hmm")
    parser.add_argument("--pfam-clans", type=Path, default=ROOT / "data" / "pfam" / "Pfam-A.clans.tsv")
    parser.add_argument("--pgap-hmm", type=Path, default=ROOT / "data" / "ncbi_pgap" / "hmm_PGAP.HMM.tgz")
    parser.add_argument("--pgap-metadata", type=Path, default=ROOT / "data" / "ncbi_pgap" / "hmm_PGAP.tsv")
    parser.add_argument("--no-pgap", action="store_true", help="Search Pfam only for the weak hits.")
    parser.add_argument("--cutoff-overrides", type=Path, default=ROOT / "data" / "pfam" / "cutoff_overrides.tsv")
    parser.add_argument("--out", type=Path, default=ROOT / "output" / "no_ga_hmm_hits.tsv")
    parser.add_argument("--cpus", type=int, default=0)
    return parser.parse_args()


def _text(value) -> str:
    return value.decode() if isinstance(value, bytes) else (value or "")


def _read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _products(genbank_dir: Path) -> dict[str, str]:
    products = {}
    for gbk in sorted(genbank_dir.glob("*.gbk")):
        for rec in SeqIO.parse(gbk, "genbank"):
            for ft in rec.features:
                if ft.type == "CDS":
                    for pid in ft.qualifiers.get("protein_id", []):
                        products[pid] = ft.qualifiers.get("product", [""])[0]
    return products


def collect_genes(segments_dir: Path) -> list[dict]:
    """Genes whose segment has no model at GA, with their group context."""
    genes = []
    without_pgap = []
    for bb_dir in sorted(p for p in segments_dir.iterdir() if (p / "groups.tsv").exists()):
        if not (bb_dir / "domains.tsv").exists():
            print(f"  skipping {bb_dir.name}: no domains.tsv (step 2 not run)", file=sys.stderr)
            continue
        params_json = bb_dir / "pfam_params.json"
        if not params_json.exists() or not json.loads(params_json.read_text()).get("pgap_hmm_path"):
            without_pgap.append(bb_dir.name)
        seg_id = {(r["sc"], r["protein_id"], r["start"], r["end"]): r["segment_id"]
                  for r in _read_tsv(bb_dir / "segments.tsv")}
        ga_rows = [r for r in _read_tsv(bb_dir / "domains.tsv")
                   if r.get("passes_ga", "True") == "True" and r["segment_id"]]
        ga_segments = {r["segment_id"] for r in ga_rows}
        group_rows = defaultdict(list)
        for r in _read_tsv(bb_dir / "groups.tsv"):
            group_rows[r["group_id"]].append(r)

        for group_id, rows in group_rows.items():
            members = {r["sc"] for r in rows}
            group_segments = {seg_id.get((r["sc"], r["protein_id"], r["start"], r["end"])) for r in rows}
            ga_scs_by_model = defaultdict(set)
            for d in ga_rows:
                if d["straddling"] == "False" and d["segment_id"] in group_segments:
                    ga_scs_by_model[f"{d['pfam_name']} [{d.get('source', 'Pfam')}]"].add(d["sc"])
            ga_summary = ", ".join(f"{name} {len(scs)}/{len(members)}"
                                   for name, scs in sorted(ga_scs_by_model.items(), key=lambda kv: (-len(kv[1]), kv[0])))
            has_signature = any(len(scs) == len(members) for scs in ga_scs_by_model.values())
            for r in rows:
                sid = seg_id.get((r["sc"], r["protein_id"], r["start"], r["end"]))
                if sid in ga_segments:
                    continue
                genes.append({
                    "bb": bb_dir.name, "group_id": group_id, "n_scs": r["n_scs"],
                    "group_has_ga_signature": has_signature, "group_ga_models": ga_summary or "-",
                    "sc": r["sc"], "protein_id": r["protein_id"], "gene_name": r["gene_name"],
                    "start": int(r["start"]), "end": int(r["end"]),
                })
    if without_pgap:
        print(f"WARNING: step 2 ran without PGAP for {len(without_pgap)} building block(s) "
              f"(e.g. {without_pgap[0]}), so genes with only PGAP hits are listed for those too. "
              "Rerun step 2 to fix.", file=sys.stderr)
    return genes


def _sequences(genes: list[dict], fasta_root: Path) -> list:
    alphabet = Alphabet.amino()
    sequences, seen = [], set()
    for bb in sorted({g["bb"] for g in genes}):
        proteins = {(p.sc_id, p.protein_id or p.local_key): p for p in load_group_proteins(fasta_root / bb)}
        for g in genes:
            key = f"{g['sc']}|{g['protein_id']}"
            if g["bb"] == bb and key not in seen:
                seen.add(key)
                p = proteins[(g["sc"], g["protein_id"])]
                sequences.append(TextSequence(name=key.encode(), sequence=p.sequence).digitize(alphabet))
    return sequences


def relaxed_search(sequences: list, hmms, source_of, n_models: int, cpus: int, hits_by_key: dict) -> None:
    """Hits without cutoffs from one library, added to hits_by_key (key -> list of hit dicts)."""
    for hits in pyhmmer.hmmsearch(hmms, sequences, cpus=cpus, E=1e6, domE=1e6, Z=n_models, domZ=n_models):
        q = hits.query
        ga_seq, ga_dom = q.cutoffs.gathering if q.cutoffs.gathering_available() else (float("nan"),) * 2
        source, family_type = None, None
        for hit in hits:
            for dom in hit.domains:
                if dom.i_evalue > MAX_I_EVALUE:
                    continue
                if source is None:
                    source, family_type = source_of(q)
                ali = dom.alignment
                hits_by_key[_text(hit.name)].append({
                    "acc": _text(q.accession).split(".")[0], "name": _text(q.name),
                    "source": source, "family_type": family_type,
                    "env_from": dom.env_from, "env_to": dom.env_to,
                    "dom_score": dom.score, "seq_score": hit.score, "i_evalue": dom.i_evalue,
                    "ga_seq": ga_seq, "ga_dom": ga_dom,
                    "hmm_cov": (ali.hmm_to - ali.hmm_from + 1) / q.M,
                })


def _frac_ga(hit: dict) -> float:
    """Domain bitscore as a fraction of the model's domain GA; -inf if the model has no GA."""
    ga = hit["ga_dom"]
    return hit["dom_score"] / ga if ga and not math.isnan(ga) and ga > 0 else float("-inf")


def main() -> int:
    args = parse_args()
    overrides = load_cutoff_overrides(args.cutoff_overrides) if args.cutoff_overrides.exists() else {}
    pfam_names = load_pfam_names(args.pfam_clans)

    genes = collect_genes(args.segments_dir)
    print(f"{len(genes)} genes without a GA-level hit from any library; searching them again without cutoffs ...")
    sequences = _sequences(genes, args.fasta_root)
    hits_by_key: dict[str, list[dict]] = defaultdict(list)
    relaxed_search(sequences, iter_hmm_file(args.pfam_hmm), lambda q: ("Pfam", ""), len(pfam_names),
                   args.cpus, hits_by_key)
    if not args.no_pgap:
        pgap = PgapLibrary.load(args.pgap_hmm, args.pgap_metadata, pfam_names)

        def pgap_source(q):
            m = pgap.model(_text(q.accession))
            return m.source, m.family_type

        relaxed_search(sequences, pgap.iter_hmms(), pgap_source, len(pgap.models), args.cpus, hits_by_key)

    products = {}
    for bb in sorted({g["bb"] for g in genes}):
        products.update(_products(args.genbank_root / bb))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_status = defaultdict(int)
    with args.out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "bb", "group_id", "n_scs", "group_has_ga_signature", "group_ga_models",
            "sc", "protein_id", "gene_name", "product", "segment", "segment_len", "status", "override_hit",
            "best_weak_source", "best_weak_acc", "best_weak_name", "best_weak_family_type", "best_weak_dom_score",
            "best_weak_ga_dom", "best_weak_frac_ga", "best_weak_seq_score", "best_weak_i_evalue", "best_weak_hmm_cov", "other_weak_hits",
        ])
        for g in genes:
            seg_from, seg_to = g["start"], g["end"]  # 1-based inclusive, as are pyhmmer env coords
            hits = []
            for h in hits_by_key.get(f"{g['sc']}|{g['protein_id']}", []):
                overlap = min(seg_to, h["env_to"]) - max(seg_from, h["env_from"]) + 1
                if overlap / (h["env_to"] - h["env_from"] + 1) >= MIN_OVERLAP:
                    hits.append(h)
            hits.sort(key=lambda h: (-_frac_ga(h), h["i_evalue"]))
            recovered = [h for h in hits if h["acc"] in overrides
                         and h["dom_score"] >= overrides[h["acc"]] and h["seq_score"] >= overrides[h["acc"]]]
            status = "recovered_by_override" if recovered else ("weak_hit_only" if hits else "no_hit")
            n_status[status] += 1
            best = hits[0] if hits else None
            o = recovered[0] if recovered else None
            w.writerow([
                g["bb"], g["group_id"], g["n_scs"], g["group_has_ga_signature"], g["group_ga_models"],
                g["sc"], g["protein_id"], g["gene_name"], products.get(g["protein_id"], ""),
                f"{seg_from}-{seg_to}", seg_to - seg_from + 1, status,
                f"{o['name']} [{o['source']}] {o['dom_score']:.1f}/{overrides[o['acc']]:g}" if o else "-",
                *([best["source"], best["acc"], best["name"], best["family_type"] or "-", f"{best['dom_score']:.1f}",
                   f"{best['ga_dom']:.1f}", f"{_frac_ga(best):.2f}", f"{best['seq_score']:.1f}", f"{best['i_evalue']:.2e}",
                   f"{best['hmm_cov']:.2f}"] if best else ["-"] * 10),
                "; ".join(f"{h['name']} [{h['source']}] {h['dom_score']:.1f}/{h['ga_dom']:.1f}" for h in hits[1:4]) or "-",
            ])

    print(f"Wrote {args.out}: " + ", ".join(f"{n} {s}" for s, n in sorted(n_status.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
