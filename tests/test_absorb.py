"""merge_uninformative_breakpoints must merge the two groups, not move a protein out of its group.

Reproduces BB0003 (L-dihydrostreptose): strD and inkA are homologous over 59-355 / 74-369 and
cluster together there, while each one's short N-terminal piece is a singleton. Collapsing the
uninformative breakpoints must leave both full-length proteins in one shared group.
"""

from rule_suggester.segments.absorb import merge_uninformative_breakpoints
from rule_suggester.segments.io_fasta import Protein
from rule_suggester.segments.segmentation import Segment


def _protein(sc, pid, gene, length):
    return Protein(key=f"{sc}|{pid}", sc_id=sc, acc="", protein_id=pid, gene_name=gene, locus_tag="",
                   local_key=pid, order=0, sequence="M" * length)


def _seg(p, start, end):
    return Segment(id=f"{p.sc_id}|{p.protein_id}|{start + 1}-{end}", key=p.key, sc_id=p.sc_id,
                   protein_id=p.protein_id, gene_name=p.gene_name, locus_tag="", start=start, end=end)


def test_fragment_singletons_fold_into_shared_group():
    strd, inka = _protein("SC0003", "BAG22759.1", "strD", 355), _protein("SC0070", "ABD59208.1", "inkA", 369)
    s1, s2 = _seg(strd, 0, 58), _seg(strd, 58, 355)
    i1, i2 = _seg(inka, 0, 73), _seg(inka, 73, 369)
    cluster_of = {s1.id: "A", i1.id: "C", s2.id: "X", i2.id: "X"}

    segments, cluster_of, merges = merge_uninformative_breakpoints([strd, inka], [s1, s2, i1, i2], cluster_of)

    assert len(merges) == 2
    assert sorted((s.gene_name, s.start, s.end) for s in segments) == [("inkA", 0, 369), ("strD", 0, 355)]
    assert len({cluster_of[s.id] for s in segments}) == 1


def test_genuine_fusion_is_left_alone():
    """A fused A+B gene in one subcluster, separate A and B in two others: groups stay apart."""
    fused = _protein("SC1", "AB", "ab", 400)
    a2, b2 = _protein("SC2", "A2", "a", 200), _protein("SC2", "B2", "b", 200)
    a3, b3 = _protein("SC3", "A3", "a", 200), _protein("SC3", "B3", "b", 200)
    f1, f2 = _seg(fused, 0, 200), _seg(fused, 200, 400)
    segs = [f1, f2, _seg(a2, 0, 200), _seg(b2, 0, 200), _seg(a3, 0, 200), _seg(b3, 0, 200)]
    cluster_of = {f1.id: "A", f2.id: "B", segs[2].id: "A", segs[3].id: "B", segs[4].id: "A", segs[5].id: "B"}

    segments, cluster_of, merges = merge_uninformative_breakpoints([fused, a2, b2, a3, b3], segs, cluster_of)

    assert merges == []
    assert len({cluster_of[s.id] for s in segments}) == 2
