"""Per-Pfam cutoff overrides replacing GA in the Pfam scan."""

import pytest
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import Background, Builder

from rule_suggester.pfam.cutoffs import load_cutoff_overrides
from rule_suggester.pfam.params import PfamParams
from rule_suggester.pfam.scan import scan_pfam_domains
from rule_suggester.segments.io_fasta import Protein

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAACRESVLAWEHPRGPGCERHHHHHH"


def _protein(seq: str = SEQ) -> Protein:
    return Protein(key="SC0001|P1", sc_id="SC0001", acc="", protein_id="P1", gene_name="g1",
                   locus_tag="", local_key="P1", order=0, sequence=seq)


def _write_hmm(path, ga: float) -> float:
    """One-sequence HMM with accession PF99999.1 and GA `ga`; returns the self-hit bitscore."""
    alphabet = Alphabet.amino()
    seq = TextSequence(name=b"TestDom", sequence=SEQ).digitize(alphabet)
    hmm, _, _ = Builder(alphabet).build(seq, Background(alphabet))
    hmm.accession = b"PF99999.1"
    hmm.cutoffs.gathering = (ga, ga)
    with open(path, "wb") as f:
        hmm.write(f)
    unfiltered = scan_pfam_domains([_protein()], PfamParams(hmm_path=path, cutoff_overrides={"PF99999": -1000.0}))
    return unfiltered[0].bitscore


def test_load_cutoff_overrides(tmp_path):
    f = tmp_path / "o.tsv"
    f.write_text("# comment\npfam_acc\tcutoff\tnote\nPF00550.32\t17.5\tPP-binding\nPF00702\t18\t\n")
    assert load_cutoff_overrides(f) == {"PF00550": 17.5, "PF00702": 18.0}
    f.write_text("pfam_acc\tcutoff\nPF00550\t17\nPF00550.1\t18\n")
    with pytest.raises(ValueError, match="more than once"):
        load_cutoff_overrides(f)


def test_override_recovers_hit_below_ga(tmp_path):
    hmm_path = tmp_path / "test.hmm"
    score = _write_hmm(hmm_path, ga=10000.0)  # GA far above any real score

    assert scan_pfam_domains([_protein()], PfamParams(hmm_path=hmm_path)) == []

    hits = scan_pfam_domains([_protein()], PfamParams(hmm_path=hmm_path, cutoff_overrides={"PF99999": score - 1}))
    assert len(hits) == 1 and hits[0].pfam_acc == "PF99999" and hits[0].passes_ga is False

    assert scan_pfam_domains([_protein()], PfamParams(hmm_path=hmm_path, cutoff_overrides={"PF99999": score + 1})) == []


def test_hit_passing_ga_is_flagged(tmp_path):
    hmm_path = tmp_path / "test.hmm"
    _write_hmm(hmm_path, ga=1.0)
    hits = scan_pfam_domains([_protein()], PfamParams(hmm_path=hmm_path))
    assert len(hits) == 1 and hits[0].passes_ga is True


def test_unknown_override_accession_raises(tmp_path):
    hmm_path = tmp_path / "test.hmm"
    _write_hmm(hmm_path, ga=1.0)
    with pytest.raises(ValueError, match="PF00001"):
        scan_pfam_domains([_protein()], PfamParams(hmm_path=hmm_path, cutoff_overrides={"PF00001": 5.0}))
