"""Hand-written rule changes (data/manual_rules.tsv) for building blocks the automatic CONDITIONS
can't express, e.g. two alternative routes to the same product, or rules to leave out.

Tab-separated, header `rule_id  field  value  note`, `#` lines are comments. field is one of:
  CONDITIONS  replaces the computed CONDITIONS (antiSMASH rule syntax, one line)
  EXTENDERS   replaces the computed EXTENDERS; value "-" removes them
  CUTOFF      replaces the estimated CUTOFF (kb)
  SUPERIORS   comma-separated rule ids this rule is inferior to; antiSMASH discards this rule's hit
              unless it credits a gene the superior's hit does not. Use it when another rule's
              CONDITIONS are a strict superset of this one's, so this rule can never fire alone
              (BB0091 hydroxymalonyl is BB0020 methoxymalonyl without Methyltransf_3).
  DROP        leaves the rule out of scripts/combine_rules.py (value ignored; say why in note)

Writing a `not X` exclusion: justify it by what this building block's reference genes *are* -- a
discrete ACP, a standalone epimerase, an adenylation domain that isn't part of an assembly line --
and say so in the note. Not by which false matches a particular test run happened to show: that
sample only contains the collisions it contains, while the rule runs against every cluster. Check
what the reference genes carry in output/segments/<group>/domains.tsv first; a false match that no
reference-gene property can exclude needs a more specific model, not a different term.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

FIELDS = ("CONDITIONS", "EXTENDERS", "CUTOFF", "SUPERIORS", "DROP")


@dataclass(frozen=True)
class ManualRule:
    rule_id: str
    conditions: str | None = None
    extenders: str | None = None  # "-" = none
    cutoff: int | None = None
    superiors: tuple[str, ...] = ()
    drop: bool = False
    notes: tuple[str, ...] = ()


def load_manual_rules(path: Path) -> dict[str, ManualRule]:
    entries: dict[str, dict] = {}
    with Path(path).open(newline="") as f:
        rows = csv.DictReader((line for line in f if line.strip() and not line.startswith("#")), delimiter="\t")
        if rows.fieldnames is None or not {"rule_id", "field", "value"} <= set(rows.fieldnames):
            raise ValueError(f"{path}: needs a header with at least 'rule_id', 'field' and 'value' columns")
        for row in rows:
            rule_id, field = row["rule_id"].strip(), row["field"].strip().upper()
            value, note = (row.get("value") or "").strip(), (row.get("note") or "").strip()
            if field not in FIELDS:
                raise ValueError(f"{path}: {rule_id}: unknown field {field!r} (expected one of {', '.join(FIELDS)})")
            e = entries.setdefault(rule_id, {"notes": []})
            key = field.lower()
            if key in e:
                raise ValueError(f"{path}: {rule_id}: {field} given more than once")
            if field == "CUTOFF":
                try:
                    e[key] = int(value)
                except ValueError:
                    raise ValueError(f"{path}: {rule_id}: CUTOFF must be a whole number of kb, got {value!r}") from None
            elif field == "SUPERIORS":
                ids = [v.strip() for v in value.replace(",", " ").split()]
                if not ids:
                    raise ValueError(f"{path}: {rule_id}: SUPERIORS needs at least one rule id")
                if len(set(ids)) != len(ids):
                    raise ValueError(f"{path}: {rule_id}: SUPERIORS has duplicate ids")
                if rule_id in ids:
                    raise ValueError(f"{path}: {rule_id}: SUPERIORS lists the rule itself")
                e[key] = tuple(ids)
            elif field == "DROP":
                e[key] = True
            else:
                if not value:
                    raise ValueError(f"{path}: {rule_id}: {field} needs a value")
                e[key] = value
            if note:
                e["notes"].append(f"{field}: {note}")
    return {
        rule_id: ManualRule(rule_id=rule_id, conditions=e.get("conditions"), extenders=e.get("extenders"),
                            cutoff=e.get("cutoff"), superiors=e.get("superiors", ()),
                            drop=e.get("drop", False), notes=tuple(e["notes"]))
        for rule_id, e in entries.items()
    }
