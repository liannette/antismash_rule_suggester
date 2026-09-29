#!/usr/bin/env python3
"""One-off diagnostic for the Dpg product group (docs/segment_grouping_spec.md, 'Real test cases'):
a 9x9 pairwise normalized-bitscore matrix among the crotonase-fold proteins DpgB/DpgC/DpgD across
all three member BGCs, to verify DIAMOND resolves them into three distinct homology blocks despite
sharing the same overall fold (rather than one large crotonase-fold blob).

Run directly against the whole proteins (not segments) since the question here is protein-level
family structure, and normalized the same way edge weights are: bitscore / min(selfbits).
"""

import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rule_suggester.headers import parse_header
from rule_suggester.segments.io_fasta import load_group_proteins

DIAMOND_BIN = Path("/Users/alien/miniconda3/envs/rule_suggester/bin/diamond")
GROUP_DIR = Path("data/fasta/SCG0031_3-5-dihydroxyphenylglycine")
OUT_PATH = Path("output/segments/SCG0031_3-5-dihydroxyphenylglycine/dpg_crotonase_matrix.tsv")

# DpgB/DpgC/DpgD orthologs in each of the three member BGCs (report.md G02/G03-G06/G07-G08).
CROTONASE_PROTEIN_IDS = {
    "CAC48379.1", "CAG15030.1", "AEI58889.1",  # DpgB
    "CAC48380.1", "CAG15031.1", "AEI58890.1",  # DpgC
    "CAC48381.1", "CAG15032.1", "AEI58891.1",  # DpgD
}
LABELS = {
    "CAC48379.1": "BGC0000311|dpgB", "CAG15030.1": "BGC0000441|dpgB(ORF18)", "AEI58889.1": "BGC0000455|vcm24",
    "CAC48380.1": "BGC0000311|dpgC", "CAG15031.1": "BGC0000441|dpgC(ORF19)", "AEI58890.1": "BGC0000455|vcm25",
    "CAC48381.1": "BGC0000311|dpgD", "CAG15032.1": "BGC0000441|dpgD(ORF20)", "AEI58891.1": "BGC0000455|vcm26",
}


def main() -> None:
    proteins = [p for p in load_group_proteins(GROUP_DIR) if p.protein_id in CROTONASE_PROTEIN_IDS]
    missing = CROTONASE_PROTEIN_IDS - {p.protein_id for p in proteins}
    if missing:
        raise SystemExit(f"missing expected protein IDs: {sorted(missing)}")

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        fasta_path = tmp_path / "crotonase.faa"
        with fasta_path.open("w") as f:
            for p in proteins:
                f.write(f">{p.header}\n{p.sequence}\n")

        db = tmp_path / "db"
        subprocess.run(
            [str(DIAMOND_BIN), "makedb", "--in", str(fasta_path), "-d", str(db)],
            check=True, capture_output=True,
        )
        hits_path = tmp_path / "hits.tsv"
        subprocess.run(
            [
                str(DIAMOND_BIN), "blastp", "-q", str(fasta_path), "-d", str(db), "-o", str(hits_path),
                "--more-sensitive", "--max-target-seqs", "0", "--max-hsps", "0", "-e", "1",
                "--outfmt", "6", "qseqid", "sseqid", "bitscore",
            ],
            check=True, capture_output=True,
        )

        best: dict[tuple[str, str], float] = {}
        selfbits: dict[str, float] = {}
        with hits_path.open() as f:
            for line in f:
                qseqid, sseqid, bitscore = line.strip().split("\t")
                q_pid = parse_header(qseqid).protein_id
                s_pid = parse_header(sseqid).protein_id
                bits = float(bitscore)
                if q_pid == s_pid:
                    selfbits[q_pid] = max(selfbits.get(q_pid, 0.0), bits)
                else:
                    key = (q_pid, s_pid)
                    best[key] = max(best.get(key, 0.0), bits)

    ordered_ids = sorted(CROTONASE_PROTEIN_IDS, key=lambda pid: (LABELS[pid].split("|")[1][:4], LABELS[pid]))

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w") as f:
        f.write("\t" + "\t".join(LABELS[pid] for pid in ordered_ids) + "\n")
        for a in ordered_ids:
            row = [LABELS[a]]
            for b in ordered_ids:
                if a == b:
                    row.append("1.000")
                else:
                    bits = best.get((a, b), best.get((b, a), 0.0))
                    denom = min(selfbits.get(a, 0.0), selfbits.get(b, 0.0))
                    row.append(f"{bits / denom:.3f}" if denom else "0.000")
            f.write("\t".join(row) + "\n")

    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
