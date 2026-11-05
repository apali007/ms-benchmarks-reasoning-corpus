import argparse
import json 
from collections import defaultdict

from llm_reasoning_corpus import spectrum_stats, add_stats, print_stats

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Print statistics for a generated reasoning corpus")
    parser.add_argument("corpus_path")
    args = parser.parse_args()

    with open(args.corpus_path) as f:
        records = json.load(f)

    stats = defaultdict(lambda: defaultdict(list))
    for record in records:
        add_stats(stats, record['dataset'], spectrum_stats(record['peptide'], record))

    print(f"Loaded {len(records)} records from {args.corpus_path}")
    print_stats(stats)