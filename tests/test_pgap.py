"""PGAP library: rule names, archive streaming, scanning, and PGAP-replaces-overlapping-Pfam selection."""

import io
import tarfile

from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import Background, Builder

from rule_suggester.pfam.libraries import PgapLibrary, iter_hmm_file, pgap_rule_name
from rule_suggester.pfam.mapping import DomainAssignment, compute_group_domains
from rule_suggester.pfam.params import PfamParams
from rule_suggester.pfam.scan import Domain, scan_pfam_domains
from rule_suggester.segments.io_fasta import Protein

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAACRESVLAWEHPRGPGCERHHHHHH"
METADATA_HEADER = "#ncbi_accession\tsource_identifier\tlabel\tsequence_cutoff\tdomain_cutoff\thmm_length\tfamily_type\tproduct_name\tsource\n"


def _hmm(name: str, accession: str, seq: str = SEQ):
    alphabet = Alphabet.amino()
    hmm, _, _ = Builder(alphabet).build(TextSequence(name=name.encode(), sequence=seq).digitize(alphabet),
                                        Background(alphabet))
    hmm.accession = accession
    hmm.cutoffs.gathering = (1.0, 1.0)
    return hmm


def _protein() -> Protein:
    return Protein(key="SC0001|P1", sc_id="SC0001", acc="", protein_id="P1", gene_name="g1",
                   locus_tag="", local_key="P1", order=0, sequence=SEQ)


def _pgap_files(tmp_path, labels: dict[str, str]):
    """Archive of one-model files (like hmm_PGAP.HMM.tgz) + metadata; labels: accession -> label."""
    archive = tmp_path / "hmm_PGAP.HMM.tgz"
    with tarfile.open(archive, "w:gz") as tf:
        for acc, label in labels.items():
            buf = io.BytesIO()
            _hmm(label, acc).write(buf)
            info = tarfile.TarInfo(f"hmm_PGAP/{acc}.HMM")
            info.size = len(buf.getvalue())
            tf.addfile(info, io.BytesIO(buf.getvalue()))
    meta = tmp_path / "hmm_PGAP.tsv"
    meta.write_text(METADATA_HEADER + "".join(
        f"{acc}\t\t{label}\t1\t1\t100\tequivalog\tproduct of {label}\tJCVI\n" for acc, label in labels.items()))
    return archive, meta


def test_pgap_rule_name():
    assert pgap_rule_name("PRK06545.2-4", "NF005112.0", set()) == "PRK06545_2-4"
    assert pgap_rule_name("AAA", "TIGR02928.1", {"AAA"}) == "AAA_TIGR02928"
    assert pgap_rule_name("FkbH", "TIGR01686.1", {"AAA"}) == "FkbH"


def test_archive_streaming_and_renaming(tmp_path):
    archive, meta = _pgap_files(tmp_path, {"TIGR00001.1": "AAA", "NF000002.1": "PRK1.2-3"})
    assert sorted(h.accession for h in iter_hmm_file(archive)) == ["NF000002.1", "TIGR00001.1"]
    lib = PgapLibrary.load(archive, meta, pfam_names={"AAA"})
    assert sorted(h.name for h in lib.iter_hmms()) == ["AAA_TIGR00001", "PRK1_2-3"]


def test_scan_pfam_and_pgap(tmp_path):
    pfam = tmp_path / "pfam.hmm"
    with pfam.open("wb") as f:
        _hmm("PfamDom", "PF99999.1").write(f)
    archive, meta = _pgap_files(tmp_path, {"TIGR00001.1": "TigrDom"})
    lib = PgapLibrary.load(archive, meta, pfam_names={"PfamDom"})
    hits = scan_pfam_domains([_protein()], PfamParams(hmm_path=pfam), lib)
    assert {(d.source, d.pfam_acc, d.pfam_name, d.family_type) for d in hits} == {
        ("Pfam", "PF99999", "PfamDom", ""), ("TIGRFAM", "TIGR00001", "TigrDom", "equivalog")}


def _assignment(sc: str, acc: str, name: str, env: tuple[int, int], source: str = "Pfam",
                family_type: str = "", i_evalue: float = 1e-10) -> DomainAssignment:
    d = Domain(sc_id=sc, protein_id=f"{sc}_p", gene_name="g", pfam_acc=acc, pfam_name=name,
               env_from=env[0], env_to=env[1], bitscore=50.0, i_evalue=i_evalue, source=source,
               family_type=family_type)
    return DomainAssignment(domain=d, segment_id=f"{sc}|seg", frac_domain_in_segment=1.0,
                            frac_segment_covered=0.5, straddling=False)


def test_pgap_replaces_overlapping_pfam_only():
    assignments = []
    for sc in ("SC1", "SC2"):
        assignments += [
            _assignment(sc, "PF00001", "ECH_1", (10, 200)),          # covered by the PGAP model
            _assignment(sc, "PF00002", "Other", (400, 500)),         # elsewhere on the protein
            _assignment(sc, "NF00001", "DpgB", (0, 300), "NCBIFAM", "equivalog", 1e-80),
            _assignment(sc, "TIGR0001", "Crotonase", (5, 290), "TIGRFAM", "subfamily", 1e-90),
        ]
    rows = compute_group_domains(assignments, {"SC1|seg": "G01", "SC2|seg": "G01"}, {"G01": 2},
                                 PfamParams(hmm_path=None))
    by_name = {g.pfam_name: g for g in rows}
    assert by_name["DpgB"].is_representative  # equivalog beats the overlapping subfamily model
    assert not by_name["Crotonase"].is_representative and by_name["Crotonase"].replaced_by == "DpgB"
    assert not by_name["ECH_1"].is_representative and by_name["ECH_1"].replaced_by == "DpgB"
    assert by_name["Other"].is_representative and by_name["Other"].replaced_by == ""


def test_prk_is_never_representative_and_replaces_nothing():
    assignments = []
    for sc in ("SC1", "SC2"):
        assignments += [
            _assignment(sc, "PF00001", "adh_short_C2", (10, 200)),
            _assignment(sc, "NF00009", "PRK12826_1-2", (0, 250), "PRK", "subfamily", 1e-90),
        ]
    rows = compute_group_domains(assignments, {"SC1|seg": "G01", "SC2|seg": "G01"}, {"G01": 2},
                                 PfamParams(hmm_path=None))
    by_name = {g.pfam_name: g for g in rows}
    assert by_name["PRK12826_1-2"].is_signature and not by_name["PRK12826_1-2"].is_representative
    assert by_name["adh_short_C2"].is_representative and by_name["adh_short_C2"].replaced_by == ""


def test_pfam_only_selection_unchanged():
    assignments = [_assignment(sc, "PF00001", "A", (0, 100)) for sc in ("SC1", "SC2")]
    rows = compute_group_domains(assignments, {"SC1|seg": "G01", "SC2|seg": "G01"}, {"G01": 2},
                                 PfamParams(hmm_path=None))
    assert [(g.pfam_name, g.is_signature, g.is_representative) for g in rows] == [("A", True, True)]
