# antismash_rule_suggester

## Setup

```bash
conda create -n rule_suggester python=3.12
conda activate rule_suggester
conda install -c bioconda diamond
pip install -e ".[dev]"
```

## Inputs

- `data/fasta/<BB…>/<SC…>.fasta` -- one directory per building block, one protein FASTA per subcluster
  (headers `sc:<SC>|acc:<accession>|pid:|gid:|locustag:`)
- `data/genbanks/<BB…>/<SC…>.gbk` -- each subcluster's reference record (cropped to its genes), needed to estimate CUTOFF
- `data/pfam/Pfam-A.hmm` and `data/pfam/Pfam-A.clans.tsv`
- `data/ncbi_pgap/hmm_PGAP.HMM.tgz` and `data/ncbi_pgap/hmm_PGAP.tsv` -- NCBI PGAP HMMs (TIGRFAM, NCBIFAM, PRK); the archive is read as is, no need to unpack it
- `data/smcogs/smcogs.hmm` -- antiSMASH's smCOG profiles (copy of `antismash/detection/genefunctions/tools/data/smcogs.hmm`), and `data/smcogs/smcog_profiles.tsv`, the curated list of profiles to use
- `data/building_blocks.json` -- curated building-block names, SMILES and chemical class

All three come from `../subcluster_collection/output/` (`build_jsons.py`, `export_fastas.py`, `export_genbanks.py`).

## Workflow

Run everything from this directory. Each group's outputs go in `output/segments/<BB…>/`; the RULE id is the `BB####` prefix.

### 1. Segment grouping (DIAMOND)

`--fasta-dir` accepts either a single product-group directory or a parent directory containing several of them:

```bash
# one group (e.g. data/fasta/BB0007_3-amino-4-7-dihydroxycoumarin)
python scripts/segment_groups.py --fasta-dir data/fasta/BB0007_3-amino-4-7-dihydroxycoumarin --out-dir output/segments

# all groups under a parent directory (e.g. data/fasta/BB0002_streptidine-6-phosphate, BB0003_L-dihydrostreptose, ...)
python scripts/segment_groups.py --fasta-dir data/fasta --out-dir output/segments
```

Writes `segments.tsv`, `groups.tsv`, `presence.tsv` and `report.md`.

### 2. Pfam/PGAP domain mapping (hmmsearch)

```bash
python scripts/pfam_domains.py --fasta-dir data/fasta --segments-dir output/segments
```

Writes `domains.tsv`, `group_domains.tsv` and `pfam_report.md`. `--fasta-dir` works the same way as in step 1.

Pfam, the NCBI PGAP library (`--no-pgap` to skip it) and the smCOG profiles listed in `data/smcogs/smcog_profiles.tsv` (columns `smcog_id`, `cutoff`, `note`; smCOGs have no GA, so the curated bitscore cutoff is used instead -- pick it from a clear score gap between target and other genes, since antiSMASH's own smCOG calls only keep each gene's best profile) are searched. The `source` column says which library a hit comes from (Pfam, TIGRFAM, NCBIFAM, PRK, smCOG), and `family_type` gives the PGAP family type (equivalog = one specific function). Representatives per group: curated smCOG signature profiles first (named `SMCOG1147` etc. in rules), then TIGRFAM and NCBIFAM signature models (most specific family type, then best E-value; overlapping ones compete); PRK models are listed but never representative, since they're full-length models with strict cutoffs that would make rules miss diverged homologs. A Pfam signature domain that an smCOG/TIGRFAM/NCBIFAM representative covers (>= 50% of its envelope) is replaced by it (`replaced_by`); the remaining Pfam signature domains keep the one-per-clan selection. PGAP models are named by their label, with characters the rule grammar doesn't allow replaced by `_` (`PRK06545.2-4` -> `PRK06545_2-4`) and the accession appended where the label is also a Pfam name (`AAA` -> `AAA_TIGR02928`).

Domains are detected with each model's gathering threshold (GA), except for Pfam or PGAP models listed in `data/pfam/cutoff_overrides.tsv` (columns `pfam_acc`, `cutoff`, `note`; another file via `--cutoff-overrides`). For those, the given bitscore replaces both GA values (sequence and domain). `domains.tsv` marks hits that only pass because of an override (`passes_ga = False`), and `pfam_report.md` lists them under "Cutoff overrides". Use the same cutoff in the antiSMASH `hmmdetails.txt` entry.

Optional, to review genes without a normal HMM hit:

```bash
python scripts/no_hmm_hits.py
```

Writes `output/no_ga_hmm_hits.tsv`: one row per gene (segment) with no model from any library (Pfam, TIGRFAM, NCBIFAM, PRK) at GA (override-only hits don't count), its group context, its best hits from a Pfam + PGAP search without cutoffs (ranked by domain bitscore / GA, closest misses first), and whether an override in `cutoff_overrides.tsv` recovers it.

### 3. Build a RULE per group

`--segments-dir` accepts a single group directory or the parent directory:

```bash
# one group (prints the full rule_conditions.md)
python scripts/build_rule.py --segments-dir output/segments/BB0050_3-5-dihydroxyphenylglycine

# all groups (prints one status line per group)
python scripts/build_rule.py --segments-dir output/segments
```

Writes `rule.txt` (the RULE block) and `rule_conditions.md` (the per-group detail and warnings behind it). GenBank files are read from `data/genbanks/<group dir name>/` by default; change the parent with `--genbank-root`, or for a single group point `--genbank-dir` straight at its GenBank directory.

- CUTOFF: for each subcluster, the smallest edge-to-edge reach that one CONDITIONS gene needs to cover all the other CONDITIONS genes (antiSMASH searches CUTOFF kb either side of an anchor CDS). The largest of these across subclusters, rounded up to the nearest 5 kb, and at least 5 kb (genes that touch in every example would otherwise give CUTOFF 0).
- NEIGHBOURHOOD is always 0.
- Complementary groups: non-core groups that never share a subcluster, together cover all of them and share a signature domain (one function split by step 1, e.g. a divergent ortholog) are merged into one core group (`G03+G04`) whose representatives are the shared signature domains.
- `data/manual_rules.tsv` (columns `rule_id`, `field`, `value`, `note`) holds hand-written changes: `CONDITIONS` replaces the computed CONDITIONS (antiSMASH syntax, one line), `EXTENDERS` replaces EXTENDERS (`-` = none), `CUTOFF` sets CUTOFF, `DROP` leaves the rule out of step 5. Use it for what the automatic rules can't express, e.g. alternative routes (`a or (b and c)`) or two paralogs that should both be required (`repeatable(2, [x])`). Manual rules are marked `# MANUAL` / `# DROPPED` in `rule.txt`.

### 4. Combined report (optional, for review)

```bash
python scripts/combined_report.py --segments-dir output/segments
```

Merges `report.md`, `pfam_report.md` and `rule.txt` into `combined_report.md`.

### 5. Combine the rules

```bash
python scripts/combine_rules.py
```

Writes `output/subcluster_rules.txt`. Groups whose `rule.txt` starts with `# WARNING` are skipped: check their `rule_conditions.md`, fix them, and re-run. Rules dropped in `data/manual_rules.tsv` are skipped too.

### 6. antiSMASH support files

```bash
python scripts/build_compound_details.py --rules-file output/subcluster_rules.txt
python scripts/extract_rule_hmms.py --rules-file output/subcluster_rules.txt --out-dir output/hmms
```

Writes `output/compound_details.txt` (one row per RULE; antiSMASH needs a row for every RULE to load the rules) and one `<domain>.hmm` per Pfam/PGAP/smCOG model the rules reference (PGAP models under their rule name), plus `cutoffs.tsv`. Each `.hmm`'s GA line is the cutoff antiSMASH should use (the same as in step 2), and `cutoffs.tsv` lists them for `hmmdetails.txt`: a cutoff override if there is one, the curated cutoff for smCOGs, otherwise the model's GA (Pfam, TIGRFAM, NCBIFAM).

### Notes

- `hmmdetails.txt` entries for the antiSMASH subclusters `data/` directory are not generated; add them by hand.
- After changing rule logic (e.g. CUTOFF), only steps 3–6 need re-running. Steps 1–2 are the slow DIAMOND/hmmsearch steps and their outputs don't change.
