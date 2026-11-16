#!/bin/bash
#SBATCH --job-name=verify_corpus
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/disk/raptor-array/apali007/msbenchmarking/slurm_logs/hardness_index_0/%x-%j.out
#SBATCH --error=/disk/raptor-array/apali007/msbenchmarking/slurm_logs/hardness_index_0/%x-%j.err

PROJECT=/disk/dragon-storage/homes/apali007/projects/ms-benchmarks
PYTHON=/disk/dragon-storage/homes/apali007/miniconda3/envs/ms-benchmarks/bin/python
CORPUS=/disk/raptor-array/apali007/msbenchmarking/hardness_index_0/reasoning_corpus/hardness_index_0.jsonl

export PYTHONUNBUFFERED=1

cd $PROJECT
date
srun $PYTHON v2/verify_corpus.py $CORPUS
date
