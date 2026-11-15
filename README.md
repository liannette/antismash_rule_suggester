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

### 2. Pfam domain mapping (hmmsearch)

```bash
python scripts/pfam_domains.py --fasta-dir data/fasta --segments-dir output/segments
```

Writes `group_domains.tsv` and `pfam_report.md`. `--fasta-dir` works the same way as in step 1.

### 3. Build a RULE per group

`--segments-dir` accepts a single group directory or the parent directory:

```bash
# one group (prints the full rule_conditions.md)
python scripts/build_rule.py --segments-dir output/segments/BB0050_3-5-dihydroxyphenylglycine

# all groups (prints one status line per group)
python scripts/build_rule.py --segments-dir output/segments
```

Writes `rule.txt` (the RULE block) and `rule_conditions.md` (the per-group detail and warnings behind it). GenBank files are read from `data/genbanks/<group dir name>/` by default; change the parent with `--genbank-root`, or for a single group point `--genbank-dir` straight at its GenBank directory.

- CUTOFF: for each subcluster, the smallest edge-to-edge reach that one CONDITIONS gene needs to cover all the other CONDITIONS genes (antiSMASH searches CUTOFF kb either side of an anchor CDS). The largest of these across subclusters, rounded up to the nearest 5 kb.
- NEIGHBOURHOOD is always 0.

### 4. Combined report (optional, for review)

```bash
python scripts/combined_report.py --segments-dir output/segments
```

Merges `report.md`, `pfam_report.md` and `rule.txt` into `combined_report.md`.

### 5. Combine the rules

```bash
python scripts/combine_rules.py
```

Writes `output/subcluster_rules.txt`. Groups whose `rule.txt` starts with `# WARNING` are skipped: check their `rule_conditions.md`, fix them, and re-run.

### 6. antiSMASH support files

```bash
python scripts/build_compound_details.py --rules-file output/subcluster_rules.txt
python scripts/extract_rule_hmms.py --rules-file output/subcluster_rules.txt --out-dir output/hmms
```

Writes `output/compound_details.txt` (one row per RULE; antiSMASH needs a row for every RULE to load the rules) and one `<domain>.hmm` per Pfam domain the rules reference.

### Notes

- `hmmdetails.txt` entries for the antiSMASH subclusters `data/` directory are not generated; add them by hand.
- After changing rule logic (e.g. CUTOFF), only steps 3–6 need re-running. Steps 1–2 are the slow DIAMOND/hmmsearch steps and their outputs don't change.
