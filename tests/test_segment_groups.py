"""Synthetic ground-truth tests for the segment-grouping pipeline (docs/segment_grouping_spec.md).

Each test builds proteins from random sequences (see synthetic.py) so homology is exact and known
in advance, then runs the real pipeline end to end and checks the written outputs.
"""

from pathlib import Path

import pandas as pd
import pytest

from rule_suggester.segments.params import Params
from rule_suggester.segments.pipeline import run_group
from synthetic import mutate, random_protein, write_sc_fasta

DIAMOND_BIN = Path("/Users/alien/miniconda3/envs/rule_suggester/bin/diamond")

pytestmark = pytest.mark.skipif(not DIAMOND_BIN.exists(), reason=f"DIAMOND not found at {DIAMOND_BIN}")


def _run(tmp_path: Path, name: str, sc_genes: dict[str, list[tuple[str, str, str]]]):
    group_dir = tmp_path / "group_in"
    group_dir.mkdir()
    for sc_id, genes in sc_genes.items():
        write_sc_fasta(group_dir, sc_id, genes)

    out_dir = tmp_path / "group_out"
    result = run_group(group_dir, out_dir, Params(), DIAMOND_BIN)
    assert result.invariant_violations == [], result.invariant_violations

    groups = pd.read_csv(out_dir / "groups.tsv", sep="\t")
    segments = pd.read_csv(out_dir / "segments.tsv", sep="\t")
    return result, groups, segments


def test_one_to_one(tmp_path):
    """3 subclusters x 3 unrelated genes -> 3 core groups of 3."""
    bases = {g: random_protein(200, seed=hash(("base", g))) for g in "ABC"}
    sc_genes = {}
    for i, sc in enumerate(["SC1", "SC2", "SC3"]):
        sc_genes[sc] = [
            (f"{gene}_{sc}", gene, mutate(bases[gene], 0.2, seed=hash((sc, gene))))
            for gene in "ABC"
        ]

    result, groups, segments = _run(tmp_path, "one_to_one", sc_genes)

    assert result.n_groups == 3
    sizes = groups.groupby("group_id")["sc"].nunique()
    assert sorted(sizes.tolist()) == [3, 3, 3]
    assert (segments["n_segments_in_protein"] == 1).all()


def test_fusion(tmp_path):
    """SC3 has gene A+B fused; groups A and B each end up with 3 members."""
    base_a = random_protein(150, seed=1)
    base_b = random_protein(150, seed=2)

    a1 = mutate(base_a, 0.2, seed=11)
    b1 = mutate(base_b, 0.2, seed=12)
    a2 = mutate(base_a, 0.2, seed=21)
    b2 = mutate(base_b, 0.2, seed=22)
    a3 = mutate(base_a, 0.2, seed=31)
    b3 = mutate(base_b, 0.2, seed=32)

    sc_genes = {
        "SC1": [("A1", "geneA", a1), ("B1", "geneB", b1)],
        "SC2": [("A2", "geneA", a2), ("B2", "geneB", b2)],
        "SC3": [("AB3", "geneAB", a3 + b3)],
    }

    result, groups, segments = _run(tmp_path, "fusion", sc_genes)

    fused = segments[segments["protein_id"] == "AB3"]
    assert len(fused) == 2, "fused protein should split into exactly 2 segments"
    junction = sorted(fused["end"].tolist())[0]
    assert abs(junction - 150) <= Params().tol

    sizes = groups.groupby("group_id")["sc"].nunique()
    assert sorted(sizes.tolist()) == [3, 3]

    fused_group_ids = set(groups[groups["sc"] == "SC3"]["group_id"])
    assert len(fused_group_ids) == 2, "the fused protein's two segments must land in two different groups"


def test_duplication(tmp_path):
    """SC2 has two close copies of gene C; they must never share a group."""
    base_c = random_protein(200, seed=100)
    c1 = mutate(base_c, 0.2, seed=101)
    c2 = mutate(base_c, 0.2, seed=102)
    c2_dup1 = mutate(c2, 0.2, seed=103)  # SC2's "ortholog" copy
    c2_dup2 = mutate(c2_dup1, 0.08, seed=104)  # recent duplicate of the copy above

    sc_genes = {
        "SC1": [("C1", "geneC", c1)],
        "SC2": [("C2a", "geneC", c2_dup1), ("C2b", "geneC", c2_dup2)],
        "SC3": [("C3", "geneC", c2)],
    }

    result, groups, segments = _run(tmp_path, "duplication", sc_genes)

    assert result.n_groups == 2
    sizes = sorted(groups.groupby("group_id")["sc"].nunique().tolist())
    assert sizes == [1, 3]

    sc2_rows = groups[groups["sc"] == "SC2"]
    assert sc2_rows["group_id"].nunique() == 2, "SC2's two copies must land in different groups"


def test_unique_gene(tmp_path):
    """A gene with no homolog anywhere becomes a full-length singleton group."""
    base_a = random_protein(180, seed=200)
    unique_d = random_protein(220, seed=201)

    sc_genes = {
        "SC1": [("A1", "geneA", mutate(base_a, 0.2, seed=211)), ("D1", "geneD", unique_d)],
        "SC2": [("A2", "geneA", mutate(base_a, 0.2, seed=212))],
        "SC3": [("A3", "geneA", mutate(base_a, 0.2, seed=213))],
    }

    result, groups, segments = _run(tmp_path, "unique_gene", sc_genes)

    d_rows = groups[groups["protein_id"] == "D1"]
    assert len(d_rows) == 1
    assert d_rows.iloc[0]["seg_len"] == len(unique_d)
    d_group_id = d_rows.iloc[0]["group_id"]
    assert (groups["group_id"] == d_group_id).sum() == 1, "D's group must be a singleton"
