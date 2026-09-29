"""Pfam clan lookup: which competing/overlapping profiles belong to the same fold family.

Used later when deciding whether two signature domains of a segment group should be ANDed
(distinct domains of one architecture) or ORed (competing profiles for the same fold) in a rule
fragment -- domains in the same clan are the latter.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Clan:
    acc: str
    name: str


def load_pfam_clans(path: Path) -> dict[str, Clan]:
    """Map Pfam accession (e.g. 'PF00001') to its Clan. Pfams with no clan assignment are omitted.

    Expects the Pfam-A.clans.tsv format (tab-separated, no header):
    pfam_acc, clan_acc, clan_name, pfam_id, description
    """
    clan_by_acc: dict[str, Clan] = {}
    with Path(path).open() as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            pfam_acc, clan_acc, clan_name = fields[0], fields[1], fields[2]
            if clan_acc:
                clan_by_acc[pfam_acc] = Clan(acc=clan_acc, name=clan_name)
    return clan_by_acc
