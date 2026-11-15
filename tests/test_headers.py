"""Header parsing and FASTA loading: sc/acc come from the header, with the file name as fallback."""

from rule_suggester.headers import make_header, parse_header
from rule_suggester.segments.io_fasta import load_group_proteins


def test_header_round_trip():
    h = parse_header(">" + make_header("SC0003", "BGC0000724", "BAG22759.1", "strD", "SGR_5930"))
    assert (h.sc_id, h.acc, h.protein_id, h.gene_name, h.locus_tag) == (
        "SC0003", "BGC0000724", "BAG22759.1", "strD", "SGR_5930"
    )


def test_blank_fields_parse_as_empty():
    h = parse_header("sc:SC0070|acc:BGC0000813|pid:ABD59208.1|gid:inkA|locustag:")
    assert h.locus_tag == ""


def test_sc_id_from_header_and_filename_fallback(tmp_path):
    (tmp_path / "SC0001.fasta").write_text(">sc:SC0001|acc:BGC0000001|pid:P1|gid:a|locustag:\nMAAA\n")
    (tmp_path / "SC0002.fasta").write_text(">pid:P2|gid:b|locustag:\nMKKK\n")  # no sc:/acc: fields
    by_pid = {p.protein_id: p for p in load_group_proteins(tmp_path)}
    assert (by_pid["P1"].sc_id, by_pid["P1"].acc) == ("SC0001", "BGC0000001")
    assert (by_pid["P2"].sc_id, by_pid["P2"].acc) == ("SC0002", "")
    assert by_pid["P1"].header == "sc:SC0001|acc:BGC0000001|pid:P1|gid:a|locustag:"
