"""Step 3: complementary-group merging, repeatable() for duplicate multi-domain terms, manual rules."""

import sys
from pathlib import Path

import pytest

from rule_suggester.rules.conditions import build_conditions, build_group_summaries
from rule_suggester.rules.manual import load_manual_rules

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import build_rule  # noqa: E402

GD_HEADER = ("group_id\tpfam_acc\tpfam_name\tclan_acc\tclan_name\tn_scs_with_domain\tn_scs_total\tprevalence\t"
             "is_signature\tany_straddling\tbest_i_evalue\tis_representative\tsource\tfamily_type\tdescription\treplaced_by\n")


def _write_group(tmp_path: Path, presence: dict[str, list[str]], domains: list[tuple]) -> Path:
    """presence: group -> subclusters; domains: (group, name, is_signature, is_representative, evalue)."""
    d = tmp_path / "BB9999_test"
    d.mkdir()
    scs = sorted({sc for v in presence.values() for sc in v})
    groups = sorted(presence)
    lines = ["sc\t" + "\t".join(groups)]
    lines += [sc + "\t" + "\t".join("1" if sc in presence[g] else "0" for g in groups) for sc in scs]
    (d / "presence.tsv").write_text("\n".join(lines) + "\n")
    rows = ["group_id\tn_scs\tsc\tprotein_id\tgene_name\tstart\tend\tseg_len"]
    for g in groups:
        for sc in presence[g]:
            rows.append(f"{g}\t{len(presence[g])}\t{sc}\t{sc}_{g}\tgene_{g}_{sc}\t1\t100\t100")
    (d / "groups.tsv").write_text("\n".join(rows) + "\n")
    (d / "group_domains.tsv").write_text(GD_HEADER + "".join(
        f"{g}\tPF0\t{name}\t{clan}\t\t1\t1\t1.0\t{sig}\tFalse\t{ev}\t{rep}\tPfam\t\t\t\n"
        for g, name, sig, rep, ev, *rest in domains for clan in [rest[0] if rest else ""]))
    return d


def test_complementary_groups_merge_into_core(tmp_path):
    d = _write_group(tmp_path, {"G01": ["SC1", "SC2", "SC3"], "G02": ["SC1", "SC2"], "G03": ["SC3"]}, [
        ("G01", "DegT", True, True, 1e-50),
        ("G02", "ADH_N", True, True, 1e-20), ("G02", "ADH_zinc_N", True, True, 1e-30),
        ("G03", "ADH_zinc_N", True, True, 1e-25),
    ])
    summaries = build_group_summaries(d / "groups.tsv", d / "presence.tsv", d / "group_domains.tsv")
    merged = [s for s in summaries if s.member_groups]
    assert len(merged) == 1 and merged[0].group_id == "G02+G03" and merged[0].is_core
    assert merged[0].representative_domains == ["ADH_zinc_N"]
    assert build_conditions(summaries).conditions == "DegT and ADH_zinc_N"


def test_merge_fallback_keeps_one_domain_per_clan(tmp_path):
    # shared Chorismate_bind/_2 (same clan) are signatures but not representatives in either group
    d = _write_group(tmp_path, {"G01": ["SC1"], "G02": ["SC2"], "G03": ["SC1", "SC2"]}, [
        ("G01", "isochor_syn", True, True, 1e-80), ("G01", "Chorismate_bind", True, False, 1e-60, "CL0969"),
        ("G01", "Chorismate_bind_2", True, False, 1e-20, "CL0969"),
        ("G02", "salicyl_syn", True, True, 1e-80), ("G02", "Chorismate_bind", True, False, 1e-55, "CL0969"),
        ("G02", "Chorismate_bind_2", True, False, 1e-25, "CL0969"),
        ("G03", "DegT", True, True, 1e-50),
    ])
    summaries = build_group_summaries(d / "groups.tsv", d / "presence.tsv", d / "group_domains.tsv")
    merged = [s for s in summaries if s.member_groups]
    assert merged[0].representative_domains == ["Chorismate_bind"]


def test_overlapping_groups_do_not_merge(tmp_path):
    d = _write_group(tmp_path, {"G01": ["SC1", "SC2", "SC3"], "G02": ["SC1", "SC2"], "G03": ["SC2", "SC3"]}, [
        ("G01", "DegT", True, True, 1e-50), ("G02", "ADH_zinc_N", True, True, 1e-30),
        ("G03", "ADH_zinc_N", True, True, 1e-25),
    ])
    summaries = build_group_summaries(d / "groups.tsv", d / "presence.tsv", d / "group_domains.tsv")
    assert not any(s.member_groups for s in summaries)


def test_duplicate_multi_domain_term_is_deduped_and_reported(tmp_path):
    d = _write_group(tmp_path, {"G01": ["SC1"], "G02": ["SC1"], "G03": ["SC1"]}, [
        ("G01", "p450", True, True, 1e-60),
        ("G02", "Dimerisation", True, True, 1e-10), ("G02", "Methyltransf_2", True, True, 1e-50),
        ("G03", "Dimerisation", True, True, 1e-12), ("G03", "Methyltransf_2", True, True, 1e-45),
    ])
    rule = build_conditions(build_group_summaries(d / "groups.tsv", d / "presence.tsv", d / "group_domains.tsv"))
    assert rule.condition_terms == ["p450", "cds(Dimerisation and Methyltransf_2)"]
    assert rule.duplicate_condition_groups == {"cds(Dimerisation and Methyltransf_2)": ["G02", "G03"]}


def test_load_manual_rules(tmp_path):
    f = tmp_path / "m.tsv"
    f.write_text("# c\nrule_id\tfield\tvalue\tnote\nBB1\tCONDITIONS\ta or (b and c)\ttwo routes\n"
                 "BB1\tEXTENDERS\t-\t\nBB2\tDROP\t\tno model\nBB3\tCUTOFF\t20\t\n")
    rules = load_manual_rules(f)
    assert rules["BB1"].conditions == "a or (b and c)" and rules["BB1"].extenders == "-"
    assert rules["BB2"].drop and rules["BB3"].cutoff == 20
    f.write_text("rule_id\tfield\tvalue\nBB1\tCONDITION\tx\n")
    with pytest.raises(ValueError, match="unknown field"):
        load_manual_rules(f)


def test_manual_conditions_and_drop_in_rule_txt(tmp_path):
    d = _write_group(tmp_path, {"G01": ["SC1"], "G02": ["SC2"], "G03": ["SC1"]}, [
        ("G01", "isochor_syn", True, True, 1e-50), ("G02", "salicyl_syn", True, True, 1e-50),
        ("G03", "CM-like", True, True, 1e-20),
    ])
    manual = tmp_path / "manual.tsv"
    manual.write_text("rule_id\tfield\tvalue\tnote\nBB9999\tCONDITIONS\tsalicyl_syn or (isochor_syn and CM-like)\troutes\n"
                      "BB9999\tEXTENDERS\t-\t\n")
    rules = load_manual_rules(manual)
    _, warnings = build_rule.build_group_rule(d, tmp_path / "no_gbk", tmp_path / "none.json", {"BB9999": "x"}, rules)
    text = (d / "rule.txt").read_text()
    assert not text.startswith("# WARNING") and warnings == []
    assert "CONDITIONS salicyl_syn or (isochor_syn and CM-like)" in text and "EXTENDERS" not in text.split("CONDITIONS")[1]

    manual.write_text("rule_id\tfield\tvalue\tnote\nBB9999\tDROP\t\tlater\n")
    build_rule.build_group_rule(d, tmp_path / "no_gbk", tmp_path / "none.json", {"BB9999": "x"}, load_manual_rules(manual))
    text = (d / "rule.txt").read_text()
    assert text.startswith("# DROPPED") and all(line.startswith("#") for line in text.strip().splitlines())
