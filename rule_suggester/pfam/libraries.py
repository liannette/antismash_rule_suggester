"""HMM libraries searched in step 2: Pfam and NCBI PGAP (TIGRFAM + NCBIFAM + PRK models).

PGAP models are renamed to their *rule name* -- the identifier used in RULE CONDITIONS/EXTENDERS
and as the <name>.hmm file name. That is the PGAP label with every character antiSMASH's rule
grammar doesn't accept replaced by "_" (PRK labels carry dots, e.g. PRK06545.2-4 -> PRK06545_2-4),
and, if that collides with a Pfam name, suffixed with the accession (AAA -> AAA_TIGR02928).
scan.py and scripts/extract_rule_hmms.py both use pgap_rule_name(), so the names always agree.
"""

from __future__ import annotations

import csv
import re
import tarfile
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from pyhmmer.plan7 import HMM, HMMFile

PGAP_SOURCES = {"JCVI": "TIGRFAM", "NCBIFAM": "NCBIFAM", "NCBI Protein Cluster (PRK)": "PRK"}
_INVALID_CHARS = re.compile(r"[^A-Za-z0-9_-]")


def _text(value: bytes | str | None) -> str:
    if value is None:
        return ""
    return value.decode() if isinstance(value, bytes) else value


def iter_hmm_file(path: Path) -> Iterator[HMM]:
    """HMMs from a plain HMMER3 file, or streamed from a .tgz/.tar.gz of one-model files (as NCBI
    distributes hmm_PGAP.HMM.tgz) without unpacking it."""
    path = Path(path)
    if path.name.endswith((".tgz", ".tar.gz")):
        with tarfile.open(path, "r:gz") as tf:
            for member in tf:
                if member.isfile() and not Path(member.name).name.startswith("."):
                    with HMMFile(tf.extractfile(member)) as hf:
                        yield from hf
    else:
        with HMMFile(path) as hf:
            yield from hf


def load_pfam_names(clans_path: Path) -> set[str]:
    """Every Pfam family name (column 4 of Pfam-A.clans.tsv, which lists all families)."""
    names = set()
    with Path(clans_path).open() as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                names.add(fields[3])
    return names


def pgap_rule_name(label: str, accession: str, pfam_names: set[str]) -> str:
    name = _INVALID_CHARS.sub("_", label)
    if name in pfam_names:
        name = f"{name}_{accession.split('.')[0]}"
    return name


@dataclass(frozen=True)
class PgapModel:
    accession: str  # with version, as in the HMM ACC line
    label: str
    source: str  # TIGRFAM / NCBIFAM / PRK
    family_type: str
    product_name: str


@dataclass(frozen=True)
class PgapLibrary:
    hmm_path: Path
    models: dict[str, PgapModel]  # accession (with version) -> model metadata
    pfam_names: set[str]

    @classmethod
    def load(cls, hmm_path: Path, metadata_path: Path, pfam_names: set[str]) -> "PgapLibrary":
        models = {}
        with Path(metadata_path).open(newline="") as f:
            header = f.readline().lstrip("#").rstrip("\n").split("\t")
            for row in csv.DictReader(f, fieldnames=header, delimiter="\t"):
                if row["source"] in PGAP_SOURCES:
                    models[row["ncbi_accession"]] = PgapModel(
                        accession=row["ncbi_accession"], label=row["label"], source=PGAP_SOURCES[row["source"]],
                        family_type=row["family_type"], product_name=row["product_name"],
                    )
        return cls(hmm_path=Path(hmm_path), models=models, pfam_names=pfam_names)

    def accessions(self) -> set[str]:
        """Accessions without version, e.g. for validating cutoff overrides."""
        return {acc.split(".")[0] for acc in self.models}

    def model(self, accession: str) -> PgapModel:
        m = self.models.get(accession)
        if m is None:
            raise ValueError(f"{accession} is in {self.hmm_path.name} but not in the PGAP metadata table")
        return m

    def iter_hmms(self) -> Iterator[HMM]:
        """PGAP HMMs, each renamed to its rule name."""
        for hmm in iter_hmm_file(self.hmm_path):
            acc = _text(hmm.accession)
            hmm.name = pgap_rule_name(self.model(acc).label, acc, self.pfam_names)
            yield hmm
