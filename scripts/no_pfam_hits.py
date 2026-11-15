#!/usr/bin/env python3
"""List every gene (segment) with no "normal" Pfam hit, i.e. no domain at Pfam's own gathering
threshold (GA), after step 2 (scripts/pfam_domains.py).

Hits that only pass because of data/pfam/cutoff_overrides.tsv (domains.tsv passes_ga = False) don't
count as normal. For each such gene, the Pfam search is rerun without cutoffs on its full-length
protein, and the hits overlapping its segment are reported with their GA, plus whether a cutoff
override now recovers the gene.

One row per gene, written to output/no_ga_pfam_hits.tsv by default. Columns:
  group_has_ga_signature  whether the group has a signature domain at GA anyway (from other genes)
  status                  recovered_by_override / weak_hit_only / no_hit
  override_hit            the best hit passing its override cutoff (name acc score/cutoff)
  best_weak_*             the highest-scoring hit overlapping the segment, whatever its score
  other_weak_hits         up to three further hits, "name score/GA"
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pyhmmer
from Bio import SeqIO
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMMFile

from rule_suggester.pfam.cutoffs import load_cutoff_overrides
from rule_suggester.segments.io_fasta import load_group_proteins

ROOT = Path(__file__).resolve().parent.parent
MIN_OVERLAP = 0.5  # share of a hit's envelope that must fall in the segment for it to count there
MAX_I_EVALUE = 1.0  # relaxed search: report domains up to this i-Evalue (Z = number of Pfam families)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--segments-dir", type=Path, default=ROOT / "output" / "segments")
    parser.add_argument("--fasta-root", type=Path, default=ROOT / "data" / "fasta")
    parser.add_argument("--genbank-root", type=Path, default=ROOT / "data" / "genbanks")
    parser.add_argument("--pfam-hmm", type=Path, default=ROOT / "data" / "pfam" / "Pfam-A.hmm")
    parser.add_argument("--pfam-clans", type=Path, default=ROOT / "data" / "pfam" / "Pfam-A.clans.tsv")
    parser.add_argument("--cutoff-overrides", type=Path, default=ROOT / "data" / "pfam" / "cutoff_overrides.tsv")
    parser.add_argument("--out", type=Path, default=ROOT / "output" / "no_ga_pfam_hits.tsv")
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
    """Genes whose segment has no domain at GA, with their group context."""
    genes = []
    for bb_dir in sorted(p for p in segments_dir.iterdir() if (p / "groups.tsv").exists()):
        if not (bb_dir / "domains.tsv").exists():
            print(f"  skipping {bb_dir.name}: no domains.tsv (step 2 not run)", file=sys.stderr)
            continue
        seg_id = {(r["sc"], r["protein_id"], r["start"], r["end"]): r["segment_id"]
                  for r in _read_tsv(bb_dir / "segments.tsv")}
        ga_segments = {r["segment_id"] for r in _read_tsv(bb_dir / "domains.tsv")
                       if r.get("passes_ga", "True") == "True" and r["segment_id"]}
        group_rows = defaultdict(list)
        for r in _read_tsv(bb_dir / "groups.tsv"):
            group_rows[r["group_id"]].append(r)
        ga_domains = defaultdict(list)  # group -> signature domains, recounted from GA-only hits
        for r in _read_tsv(bb_dir / "group_domains.tsv"):
            ga_domains[r["group_id"]].append(r)

        for group_id, rows in group_rows.items():
            members = {r["sc"] for r in rows}
            ga_scs_by_pfam = defaultdict(set)
            for d in _read_tsv(bb_dir / "domains.tsv"):
                if d.get("passes_ga", "True") == "True" and d["straddling"] == "False":
                    for r in rows:
                        if seg_id.get((r["sc"], r["protein_id"], r["start"], r["end"])) == d["segment_id"]:
                            ga_scs_by_pfam[d["pfam_name"]].add(d["sc"])
            ga_summary = ", ".join(f"{name} {len(scs)}/{len(members)}"
                                   for name, scs in sorted(ga_scs_by_pfam.items(), key=lambda kv: (-len(kv[1]), kv[0])))
            has_signature = any(len(scs) == len(members) for scs in ga_scs_by_pfam.values())
            for r in rows:
                sid = seg_id.get((r["sc"], r["protein_id"], r["start"], r["end"]))
                if sid in ga_segments:
                    continue
                genes.append({
                    "bb": bb_dir.name, "group_id": group_id, "n_scs": r["n_scs"],
                    "group_has_ga_signature": has_signature, "group_ga_domains": ga_summary or "-",
                    "sc": r["sc"], "protein_id": r["protein_id"], "gene_name": r["gene_name"],
                    "start": int(r["start"]), "end": int(r["end"]),
                })
    return genes


def relaxed_search(genes: list[dict], fasta_root: Path, pfam_hmm: Path, n_families: int, cpus: int):
    """Pfam hits without cutoffs on the genes' full-length proteins: key -> list of hit dicts."""
    alphabet = Alphabet.amino()
    sequences, seen = [], set()
    for bb in sorted({g["bb"] for g in genes}):
        proteins = {(p.sc_id, p.protein_id or p.local_key): p for p in load_group_proteins(fasta_root / bb)}
        for g in genes:
            if g["bb"] != bb:
                continue
            key = f"{g['sc']}|{g['protein_id']}"
            if key not in seen:
                seen.add(key)
                p = proteins[(g["sc"], g["protein_id"])]
                sequences.append(TextSequence(name=key.encode(), sequence=p.sequence).digitize(alphabet))

    hits_by_key = defaultdict(list)
    with HMMFile(pfam_hmm) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, sequences, cpus=cpus, E=1e6, domE=1e6, Z=n_families, domZ=n_families):
            q = hits.query
            ga_seq, ga_dom = q.cutoffs.gathering if q.cutoffs.gathering_available() else (float("nan"),) * 2
            for hit in hits:
                for dom in hit.domains:
                    if dom.i_evalue > MAX_I_EVALUE:
                        continue
                    ali = dom.alignment
                    hits_by_key[_text(hit.name)].append({
                        "acc": _text(q.accession).split(".")[0], "name": _text(q.name),
                        "env_from": dom.env_from, "env_to": dom.env_to,
                        "dom_score": dom.score, "seq_score": hit.score, "i_evalue": dom.i_evalue,
                        "ga_seq": ga_seq, "ga_dom": ga_dom,
                        "hmm_cov": (ali.hmm_to - ali.hmm_from + 1) / q.M,
                    })
    return hits_by_key


def main() -> int:
    args = parse_args()
    overrides = load_cutoff_overrides(args.cutoff_overrides) if args.cutoff_overrides.exists() else {}
    n_families = sum(1 for line in args.pfam_clans.open() if line.strip())

    genes = collect_genes(args.segments_dir)
    print(f"{len(genes)} genes without a GA-level Pfam hit; rerunning Pfam without cutoffs on them ...")
    hits_by_key = relaxed_search(genes, args.fasta_root, args.pfam_hmm, n_families, args.cpus)
    products = {}
    for bb in sorted({g["bb"] for g in genes}):
        products.update(_products(args.genbank_root / bb))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "bb", "group_id", "n_scs", "group_has_ga_signature", "group_ga_domains",
            "sc", "protein_id", "gene_name", "product", "segment", "segment_len", "status", "override_hit",
            "best_weak_acc", "best_weak_name", "best_weak_dom_score", "best_weak_ga_dom", "best_weak_seq_score",
            "best_weak_i_evalue", "best_weak_hmm_cov", "other_weak_hits",
        ])
        n_status = defaultdict(int)
        for g in genes:
            seg_from, seg_to = g["start"], g["end"]  # 1-based inclusive, as are pyhmmer env coords
            hits = []
            for h in hits_by_key.get(f"{g['sc']}|{g['protein_id']}", []):
                overlap = min(seg_to, h["env_to"]) - max(seg_from, h["env_from"]) + 1
                if overlap / (h["env_to"] - h["env_from"] + 1) >= MIN_OVERLAP:
                    hits.append(h)
            hits.sort(key=lambda h: -h["dom_score"])
            recovered = [h for h in hits if h["acc"] in overrides
                         and h["dom_score"] >= overrides[h["acc"]] and h["seq_score"] >= overrides[h["acc"]]]
            status = "recovered_by_override" if recovered else ("weak_hit_only" if hits else "no_hit")
            n_status[status] += 1
            best = hits[0] if hits else None
            o = recovered[0] if recovered else None
            w.writerow([
                g["bb"], g["group_id"], g["n_scs"], g["group_has_ga_signature"], g["group_ga_domains"],
                g["sc"], g["protein_id"], g["gene_name"], products.get(g["protein_id"], ""),
                f"{seg_from}-{seg_to}", seg_to - seg_from + 1, status,
                f"{o['name']} {o['acc']} {o['dom_score']:.1f}/{overrides[o['acc']]:g}" if o else "-",
                *( [best["acc"], best["name"], f"{best['dom_score']:.1f}", f"{best['ga_dom']:.1f}",
                    f"{best['seq_score']:.1f}", f"{best['i_evalue']:.2e}", f"{best['hmm_cov']:.2f}"]
                   if best else ["-"] * 7 ),
                "; ".join(f"{h['name']} {h['dom_score']:.1f}/{h['ga_dom']:.1f}" for h in hits[1:4]) or "-",
            ])

    print(f"Wrote {args.out}: " + ", ".join(f"{n} {s}" for s, n in sorted(n_status.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
