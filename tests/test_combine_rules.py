"""Step 5: combining rule.txt blocks, and the SUPERIORS ordering antiSMASH requires."""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import combine_rules  # noqa: E402


def _block(rule_id: str, superiors: str = "") -> str:
    sup = f"    SUPERIORS {superiors}\n" if superiors else ""
    return (f"RULE {rule_id}\n    CATEGORY Subcluster\n    DESCRIPTION {rule_id}\n"
            f"{sup}    CUTOFF 5\n    NEIGHBOURHOOD 0\n    CONDITIONS model_{rule_id}\n")


def _run(tmp_path: Path, blocks: dict[str, str], extra: tuple[str, ...] = ()) -> tuple[int, Path]:
    segments = tmp_path / "segments"
    for name, text in blocks.items():
        d = segments / name
        d.mkdir(parents=True)
        (d / "rule.txt").write_text(text)
    out = tmp_path / "rules.txt"
    argv = ["combine_rules.py", "--segments-dir", str(segments), "--out", str(out), *extra]
    old = sys.argv
    sys.argv = argv
    try:
        return combine_rules.main(), out
    finally:
        sys.argv = old


def test_superior_is_written_before_its_inferior(tmp_path):
    # ids sort BB0084 < BB0098, but BB0098 is the superior, so it has to come first
    code, out = _run(tmp_path, {"BB0084_x": _block("BB0084", "BB0098"), "BB0098_y": _block("BB0098")},
                     ("--no-comments",))
    assert code == 0
    order = re.findall(r"(?m)^RULE (\S+)", out.read_text())
    assert order == ["BB0098", "BB0084"]


def test_unrelated_rules_keep_id_order(tmp_path):
    code, out = _run(tmp_path, {f"BB000{n}_x": _block(f"BB000{n}") for n in (3, 1, 2)}, ("--no-comments",))
    assert code == 0
    assert re.findall(r"(?m)^RULE (\S+)", out.read_text()) == ["BB0001", "BB0002", "BB0003"]


def test_dangling_superior_is_an_error(tmp_path):
    # the superior is dropped, so it never reaches the combined file
    code, out = _run(tmp_path, {"BB0084_x": _block("BB0084", "BB0098"),
                                "BB0098_y": "# DROPPED (data/manual_rules.tsv): later\n"})
    assert code == 1 and not out.exists()


def test_superiors_cycle_is_an_error(tmp_path):
    code, _ = _run(tmp_path, {"BB1_a": _block("BB1", "BB2"), "BB2_b": _block("BB2", "BB1")})
    assert code == 1
