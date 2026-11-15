"""Reading identifiers and rule ids from a combined rules file."""

from rule_suggester.export.parse import extract_domain_ids, extract_rule_ids

RULES = """# Auto-generated header
RULE BB0001
    CATEGORY Subcluster
    DESCRIPTION something
    CUTOFF 10
    NEIGHBOURHOOD 0
    CONDITIONS DegT_DnrJ_EryC1
               and repeatable(2, [Methyltransf_2])
    EXTENDERS cds(PP-binding and AMP-binding_C)

# MANUAL (data/manual_rules.tsv): CONDITIONS: two routes: salicylate synthase (mbtI) or pchA + pchB
RULE BB0062
    CATEGORY Subcluster
    DESCRIPTION salicylic acid
    CUTOFF 5
    NEIGHBOURHOOD 0
    CONDITIONS salicyl_syn or (isochor_syn and CM-like)  # trailing comment with words
# RULE BB0123 (dropped, commented out)
#     CONDITIONS Methyltransf_25
"""


def test_identifiers_from_conditions_and_extenders_only():
    assert extract_domain_ids(RULES) == {
        "DegT_DnrJ_EryC1", "Methyltransf_2", "PP-binding", "AMP-binding_C", "salicyl_syn", "isochor_syn", "CM-like",
    }


def test_rule_ids_ignore_commented_rules():
    assert extract_rule_ids(RULES) == ["BB0001", "BB0062"]
