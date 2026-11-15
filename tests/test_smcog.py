"""Curated smCOG profiles: loading, renaming, cutoffs as GA, and their priority as representatives."""

import pytest

from rule_suggester.pfam.libraries import SmcogLibrary, load_smcog_profiles
from rule_suggester.pfam.mapping import compute_group_domains
from rule_suggester.pfam.params import PfamParams

from test_pgap import _assignment, _hmm


def _smcog_file(tmp_path):
    path = tmp_path / "smcogs.hmm"
    with path.open("wb") as f:
        for name in ("SMCOG1147:putative_acyl_carrier_protein", "SMCOG1119:halogenase", "SMCOG1000:ABC_transporter"):
            h = _hmm(name, "")
            h.accession = None  # smCOGs ship without accession ...
            h.cutoffs.gathering = None  # ... and without GA
            h.write(f)
    return path


def test_load_smcog_profiles(tmp_path):
    f = tmp_path / "p.tsv"
    f.write_text("# comment\nsmcog_id\tcutoff\tnote\nSMCOG1147\t70\tACP\nsmcog1119\t200\t\n")
    assert load_smcog_profiles(f) == {"SMCOG1147": 70.0, "SMCOG1119": 200.0}


def test_only_listed_profiles_renamed_with_cutoff_as_ga(tmp_path):
    lib = SmcogLibrary(_smcog_file(tmp_path), {"SMCOG1147": 70.0, "SMCOG1119": 200.0})
    hmms = {h.name: h for h in lib.iter_hmms()}
    assert set(hmms) == {"SMCOG1147", "SMCOG1119"}
    assert hmms["SMCOG1147"].accession == "SMCOG1147"
    assert hmms["SMCOG1147"].description == "putative acyl carrier protein"
    assert hmms["SMCOG1119"].cutoffs.gathering == (200.0, 200.0)


def test_unknown_profile_raises(tmp_path):
    lib = SmcogLibrary(_smcog_file(tmp_path), {"SMCOG9999": 50.0})
    with pytest.raises(ValueError, match="SMCOG9999"):
        list(lib.iter_hmms())


def test_smcog_representative_before_pgap_and_pfam():
    assignments = []
    for sc in ("SC1", "SC2"):
        assignments += [
            _assignment(sc, "PF00550", "PP-binding", (5, 70)),
            _assignment(sc, "TIGR0001", "ACP_x", (0, 80), "TIGRFAM", "equivalog", 1e-30),
            _assignment(sc, "SMCOG1147", "SMCOG1147", (0, 80), "smCOG", "", 1e-20),
        ]
    rows = compute_group_domains(assignments, {"SC1|seg": "G01", "SC2|seg": "G01"}, {"G01": 2},
                                 PfamParams(hmm_path=None))
    by_name = {g.pfam_name: g for g in rows}
    assert by_name["SMCOG1147"].is_representative
    assert by_name["ACP_x"].replaced_by == "SMCOG1147" and not by_name["ACP_x"].is_representative
    assert by_name["PP-binding"].replaced_by == "SMCOG1147"
