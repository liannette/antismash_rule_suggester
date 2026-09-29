# antismash_rule_suggester

## Setup

```bash
conda create -n rule_suggester python=3.12
conda activate rule_suggester
conda install -c bioconda diamond
pip install -e ".[dev]"
```

## Usage

`--fasta-dir` accepts either a single product-group directory or a parent directory containing several of them:

```bash
# one group (e.g. data/fasta/SCG0006_3-amino-4-7-dihydroxycoumarin)
python scripts/segment_groups.py --fasta-dir data/fasta/SCG0006_3-amino-4-7-dihydroxycoumarin --out-dir output/segments

# all groups under a parent directory (e.g. data/fasta/SCG0001_streptidine, SCG0002_L-streptose, ...)
python scripts/segment_groups.py --fasta-dir data/fasta --out-dir output/segments
```
