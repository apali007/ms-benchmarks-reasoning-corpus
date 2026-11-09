import argparse
import json 
from collections import defaultdict

from llm_reasoning_corpus import spectrum_stats, add_stats, print_stats

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Print statistics for a generated reasoning corpus")
    parser.add_argument("corpus_path")
    args = parser.parse_args()

    stats = defaultdict(lambda: defaultdict(list))
    n = 0
    with open(args.corpus_path) as f:
        for line in f:
            record = json.loads(line)
            add_stats(stats, record['dataset'], spectrum_stats(record['modified_peptide'], record))
            n += 1

    print(f"Loaded {n} records from {args.corpus_path}")
    print_stats(stats)