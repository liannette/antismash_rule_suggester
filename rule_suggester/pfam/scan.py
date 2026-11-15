"""Running Pfam hmmsearch against full-length proteins (spec Phase 2, step 1).

Runs on the original full-length proteins, never on segments -- segment boundaries are only
accurate to +/- TOL and would truncate domains that straddle them.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMM, HMMFile

from ..segments.io_fasta import Protein
from .params import PfamParams


@dataclass(frozen=True)
class Domain:
    sc_id: str
    protein_id: str
    gene_name: str
    pfam_acc: str
    pfam_name: str
    env_from: int  # 0-based half-open, matching the rest of the pipeline
    env_to: int
    bitscore: float
    i_evalue: float
    passes_ga: bool = True  # False: detected only because of a cutoff override below Pfam's GA


def _text(value: bytes | str | None) -> str:
    if value is None:
        return ""
    return value.decode() if isinstance(value, bytes) else value


def _with_overrides(hmms: Iterable[HMM], overrides: dict[str, float], original_ga: dict[str, tuple[float, float]]):
    """Yield each HMM, replacing its gathering thresholds (sequence and domain) with the override
    for its accession, if any. The original GA of overridden HMMs is recorded in original_ga.
    """
    for hmm in hmms:
        acc = _text(hmm.accession).split(".")[0]
        if acc in overrides:
            original_ga[acc] = hmm.cutoffs.gathering
            hmm.cutoffs.gathering = (overrides[acc], overrides[acc])
        yield hmm


def scan_pfam_domains(proteins: list[Protein], params: PfamParams) -> list[Domain]:
    """Run hmmsearch with Pfam's gathering-threshold cutoffs (--cut_ga) against every protein,
    except for Pfams in params.cutoff_overrides, which use that bitscore instead.
    Proteins are keyed by "<sc_id>|<protein_id>" for the search; that's decoded back into
    per-domain sc_id/protein_id in the result, no other header parsing involved.
    """
    alphabet = Alphabet.amino()
    proteins_by_key: dict[str, Protein] = {}
    sequences = []
    for p in proteins:
        key = f"{p.sc_id}|{p.protein_id or p.local_key}"
        proteins_by_key[key] = p
        seq = TextSequence(name=key.encode(), sequence=p.sequence).digitize(alphabet)
        sequences.append(seq)

    domains: list[Domain] = []
    original_ga: dict[str, tuple[float, float]] = {}
    with HMMFile(params.hmm_path) as hmm_file:
        queries = _with_overrides(hmm_file, params.cutoff_overrides, original_ga)
        for hits in pyhmmer.hmmsearch(queries, sequences, cpus=params.cpus, bit_cutoffs="gathering"):
            pfam_acc = _text(hits.query.accession).split(".")[0]
            ga_seq, ga_dom = original_ga.get(pfam_acc, (float("-inf"), float("-inf")))
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
                            sc_id=protein.sc_id,
                            protein_id=protein.protein_id or protein.local_key,
                            gene_name=protein.gene_name,
                            pfam_acc=pfam_acc,
                            pfam_name=pfam_name,
                            env_from=domain.env_from - 1,
                            env_to=domain.env_to,
                            bitscore=domain.score,
                            i_evalue=domain.i_evalue,
                            passes_ga=hit.score >= ga_seq and domain.score >= ga_dom,
                        )
                    )

    missing = sorted(set(params.cutoff_overrides) - set(original_ga))
    if missing:
        raise ValueError(f"cutoff override(s) for Pfam(s) not in {params.hmm_path.name}: {', '.join(missing)}")
    return domains
