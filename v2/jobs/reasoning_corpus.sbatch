#!/bin/bash
#SBATCH --job-name=reasoning_corpus
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/disk/raptor-array/apali007/msbenchmarking/hardness_index_0/reasoning_corpus/%x-%j.out
#SBATCH --error=/disk/raptor-array/apali007/msbenchmarking/hardness_index_0/reasoning_corpus/%x-%j.err

PROJECT=/disk/dragon-storage/homes/apali007/projects/ms-benchmarks
PYTHON=/disk/raptor-array/mshab008/msbenchmarking/.venv/bin/python
INPUT=/disk/raptor-array/mshab008/msbenchmarking/hardness_index_0/filtered_spectra/hardness_index_0.parquet
OUT=/disk/raptor-array/apali007/msbenchmarking/hardness_index_0/reasoning_corpus

export PYTHONUNBUFFERED=1

cd $PROJECT
date
srun $PYTHON v2/llm_reasoning_corpus.py \
    --input_data_path $INPUT \
    --num_examples 0 \
    --output_dir $OUT
date
