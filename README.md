# MS Benchmarks: LLM Reasoning Corpus

Generates step-by-step de novo peptide sequencing traces (b/y ion ladders)
from labeled PSMs, for training LLMs. Input is the `hardness_index_0` PSM
parquet file.

## Layout

```
v1/                               first prototype (spectrum descriptions), not maintained
v2/
├── llm_reasoning_corpus.py       main script: parquet → reasoning corpus (.jsonl) + stats
├── constants.py                  residue / PTM masses, data paths
├── corpus_stats.py               re-print stats from a generated .jsonl
├── analyse_ptms.py               inventory of PTMs in a PSM parquet
├── analyse_parquet.py            quick look at the parquet schema and rows
├── generate_spectra_reasoning_corpus.py   earlier version of the generator
├── jobs/reasoning_corpus.sbatch  Slurm job for the full corpus
└── output/                       local test runs (git-ignored)
```

## How a trace is built (`llm_reasoning_corpus.py`)

1. `parse_parquet`: read spectra, peptides and modified peptides.
2. `load_search_params`: per-dataset fragment tolerance and fixed mods from the
   metadata JSON; `load_variable_mods`: allowed variable mods from Tide params.
3. `tokenize_peptide` / `add_modified_residues`: split the modified peptide
   into residues, each with its own mass (e.g. `M[15.9949]`).
4. `match_ions`: match observed peaks to b/y ions (charges 1..z, H2O/NH3 losses).
5. `derive_peptide` / `extend_ladder`: write the y ladder, then the b ladder.
   Residues are named only from the first ion, a mass difference, or the
   precursor; PTMs are named on the first ion that carries them.
6. `spectrum_stats` / `print_stats`: ion and residue coverage, fragment
   error, per dataset and modified vs unmodified.

## Run

Test locally:
```bash
python v2/llm_reasoning_corpus.py --num_examples 100   # writes v2/output/hardness_index_0.jsonl
```

Full corpus:
```bash
sbatch v2/jobs/reasoning_corpus.sbatch
```

Re-print stats from an existing corpus:
```bash
python v2/corpus_stats.py <corpus.jsonl>
```

## Output record (one JSON object per line)

`id`, `dataset`, `input` (precursor, mods, tolerance, peak list), `output`
(the reasoning trace), `peptide`, `modified_peptide`, `y_ions` / `b_ions`
(matched peaks with theoretical m/z, charge, loss), `y_residues` /
`b_residues` (position and how each residue was identified).
