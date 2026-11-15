"""Running hmmsearch against full-length proteins (spec Phase 2, step 1): Pfam, plus the NCBI PGAP
library (TIGRFAM/NCBIFAM/PRK) and a curated set of antiSMASH smCOG profiles when given.

Runs on the original full-length proteins, never on segments -- segment boundaries are only
accurate to +/- TOL and would truncate domains that straddle them.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMM

from ..segments.io_fasta import Protein
from .libraries import PgapLibrary, SmcogLibrary, iter_hmm_file
from .params import PfamParams


@dataclass(frozen=True)
class Domain:
    sc_id: str
    protein_id: str
    gene_name: str
    pfam_acc: str  # model accession without version (PF00550, TIGR01686, NF042431)
    pfam_name: str  # rule name: Pfam NAME, or the PGAP rule name (see libraries.pgap_rule_name)
    env_from: int  # 0-based half-open, matching the rest of the pipeline
    env_to: int
    bitscore: float
    i_evalue: float
    passes_ga: bool = True  # False: detected only because of a cutoff override below the model's GA
    source: str = "Pfam"  # Pfam / TIGRFAM / NCBIFAM / PRK / smCOG
    family_type: str = ""  # PGAP family type (equivalog, subfamily, ...); empty for Pfam
    description: str = ""


@dataclass(frozen=True)
class _ModelInfo:
    source: str
    family_type: str
    description: str


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


def _search(
    hmms: Iterable[HMM],
    sequences: list,
    proteins_by_key: dict[str, Protein],
    params: PfamParams,
    info_of: Callable[[HMM], _ModelInfo],
    original_ga: dict[str, tuple[float, float]],
) -> list[Domain]:
    domains: list[Domain] = []
    queries = _with_overrides(hmms, params.cutoff_overrides, original_ga)
    for hits in pyhmmer.hmmsearch(queries, sequences, cpus=params.cpus, bit_cutoffs="gathering"):
        query = hits.query
        pfam_acc = _text(query.accession).split(".")[0]
        ga_seq, ga_dom = original_ga.get(pfam_acc, (float("-inf"), float("-inf")))
        info = None
        for hit in hits:
            if not hit.included:
                continue
            info = info or info_of(query)
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
                        pfam_name=_text(query.name),
                        env_from=domain.env_from - 1,
                        env_to=domain.env_to,
                        bitscore=domain.score,
                        i_evalue=domain.i_evalue,
                        passes_ga=hit.score >= ga_seq and domain.score >= ga_dom,
                        source=info.source,
                        family_type=info.family_type,
                        description=info.description,
                    )
                )
    return domains


def scan_pfam_domains(
    proteins: list[Protein], params: PfamParams, pgap: PgapLibrary | None = None, smcog: SmcogLibrary | None = None
) -> list[Domain]:
    """Run hmmsearch with each model's gathering-threshold cutoffs (--cut_ga) against every protein,
    except for models in params.cutoff_overrides, which use that bitscore instead. Pfam always;
    the PGAP library and the curated smCOG profiles (with their curated cutoffs as GA) too when given.
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

    original_ga: dict[str, tuple[float, float]] = {}
    domains = _search(
        iter_hmm_file(params.hmm_path), sequences, proteins_by_key, params,
        lambda hmm: _ModelInfo("Pfam", "", _text(hmm.description)), original_ga,
    )
    libraries = params.hmm_path.name
    if pgap is not None:
        def pgap_info(hmm: HMM) -> _ModelInfo:
            m = pgap.model(_text(hmm.accession))
            return _ModelInfo(m.source, m.family_type, m.product_name)

        domains += _search(pgap.iter_hmms(), sequences, proteins_by_key, params, pgap_info, original_ga)
        libraries += f" / {pgap.hmm_path.name}"
    if smcog is not None:
        domains += _search(smcog.iter_hmms(), sequences, proteins_by_key, params,
                           lambda hmm: _ModelInfo("smCOG", "", _text(hmm.description)), original_ga)
        libraries += f" / {smcog.hmm_path.name}"

    missing = sorted(set(params.cutoff_overrides) - set(original_ga))
    if missing:
        raise ValueError(f"cutoff override(s) for model(s) not in {libraries}: {', '.join(missing)}")
    return domains
