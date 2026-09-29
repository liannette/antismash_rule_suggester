"""Running Pfam hmmsearch against full-length proteins (spec Phase 2, step 1).

Runs on the original full-length proteins, never on segments -- segment boundaries are only
accurate to +/- TOL and would truncate domains that straddle them.
"""

from __future__ import annotations

from dataclasses import dataclass

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMMFile

from ..segments.io_fasta import Protein
from .params import PfamParams


@dataclass(frozen=True)
class Domain:
    bgc_id: str
    protein_id: str
    gene_name: str
    pfam_acc: str
    pfam_name: str
    env_from: int  # 0-based half-open, matching the rest of the pipeline
    env_to: int
    bitscore: float
    i_evalue: float


def _text(value: bytes | str | None) -> str:
    if value is None:
        return ""
    return value.decode() if isinstance(value, bytes) else value


def scan_pfam_domains(proteins: list[Protein], params: PfamParams) -> list[Domain]:
    """Run hmmsearch with Pfam's gathering-threshold cutoffs (--cut_ga) against every protein.
    Proteins are keyed by "<bgc_id>|<protein_id>" for the search; that's decoded back into
    per-domain bgc_id/protein_id in the result, no other header parsing involved.
    """
    alphabet = Alphabet.amino()
    proteins_by_key: dict[str, Protein] = {}
    sequences = []
    for p in proteins:
        key = f"{p.bgc_id}|{p.protein_id or p.local_key}"
        proteins_by_key[key] = p
        seq = TextSequence(name=key.encode(), sequence=p.sequence).digitize(alphabet)
        sequences.append(seq)

    domains: list[Domain] = []
    with HMMFile(params.hmm_path) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, sequences, cpus=params.cpus, bit_cutoffs="gathering"):
            pfam_acc = _text(hits.query.accession).split(".")[0]
            pfam_name = _text(hits.query.name)
            for hit in hits:
                if not hit.included:
                    continue
                protein = proteins_by_key[_text(hit.name)]
                for domain in hit.domains:
                    if not domain.included:
                        continue
                    domains.append(
                        Domain(
                            bgc_id=protein.bgc_id,
                            protein_id=protein.protein_id or protein.local_key,
                            gene_name=protein.gene_name,
                            pfam_acc=pfam_acc,
                            pfam_name=pfam_name,
                            env_from=domain.env_from - 1,
                            env_to=domain.env_to,
                            bitscore=domain.score,
                            i_evalue=domain.i_evalue,
                        )
                    )
    return domains
