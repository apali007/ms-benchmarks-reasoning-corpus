import re
import argparse
from collections import Counter

import pyarrow.parquet as pq

# One residue plus any modifications written after it, e.g M[15.9949]
TOKEN = re.compile(r'([A-Z])((?:\[[^\]]+\])*)')
DELTA = re.compile(r'\[([^\]]+)\]')

def position_name(index, length):
    if index == 0:
        return 'first'
    return 'last' if index == length -1 else 'internal'

def scan_peptide(modified_peptide):
    # Returns (residue, delta, position) for every modification in the peptide
    tokens = TOKEN.findall(modified_peptide)
    sites = []
    for index, (residue, mods) in enumerate(tokens):
        for delta in DELTA.findall(mods):
            sites.append((residue, delta, position_name(index, len(tokens))))
    
    return sites

def analyse(parquet_path):
    columns = ['dataset', 'peptide', 'modified_peptide']
    total, modified, mods_per_peptide = Counter(), Counter(), Counter()
    sites_found, unparsed = Counter(), []

    for batch in pq.ParquetFile(parquet_path).iter_batches(batch_size=100000, columns=columns):
        df = batch.to_pandas()
        for dataset, peptide, modified_peptide in zip(df['dataset'], df['peptide'], df['modified_peptide']):
            total[dataset] += 1
            if peptide == modified_peptide:
                continue

            modified[dataset] += 1

            # A format this script does not understand must be reported, not silently miscounted
            sites = scan_peptide(modified_peptide)
            if DELTA.sub('', modified_peptide) != peptide or len(sites) != len(DELTA.findall(modified_peptide)):
                unparsed.append(modified_peptide)
                continue

            mods_per_peptide[len(sites)] += 1
            for residue, delta, position in sites:
                sites_found[(dataset, delta, residue, position)] += 1

    return total, modified, mods_per_peptide, sites_found, unparsed

def print_report(total, modified, mods_per_peptide, sites_found, unparsed):
    
    print("Modified PSMs per dataset:")
    print(f"{'dataset':10s} {'total':>9s} {'modified':>9s} {'percent':>8s}")
    
    for dataset in sorted(total):
        percent = 100 * modified[dataset] / total[dataset]
        print(f"{dataset:10s} {total[dataset]:9d} {modified[dataset]:9d} {percent:7.2f}%")
    
    percent = 100 * sum(modified.values()) / sum(total.values())
    print(f"{'all':10s} {sum(total.values()):9d} {sum(modified.values()):9d} {percent:7.2f}%")

    print("\nModification sites:")
    print(f"{'dataset':10s} {'delta':>10s} {'residue':>8s} {'position':>9s} {'sites':>7s}")
    for (dataset, delta, residue, position), count in sorted(sites_found.items()):
        print(f"{dataset:10s} {delta:>10s} {residue:>8s} {position:>9s} {count:7d}")

    print("\nModifications per peptide:")
    for n_mods, count in sorted(mods_per_peptide.items()):
        print(f"{n_mods}: {count}")

    print(f"\nUnparsed peptides: {len(unparsed)}")
    for peptide in unparsed[:10]:
        print(peptide)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Report which PTMs are present in a PSM parquet file")
    parser.add_argument("--input_data_path", 
                        default="/disk/raptor-array/mshab008/msbenchmarking/hardness_index_0/filtered_spectra/hardness_index_0.parquet")

    args = parser.parse_args()

    print_report(*analyse(args.input_data_path))


    