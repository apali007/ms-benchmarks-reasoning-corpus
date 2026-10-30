import json
import re
import numpy as np
import pyarrow.parquet as pq


AA_MASSES = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}

MASS_H2O    = 18.01056
MASS_PROTON =  1.00728
MASS_C13    =  1.00335

# Compute theoretical b and y ions for a peptide sequence
def theoretical_ions(peptide, charge=1):
    b_ions, y_ions = {}, {}
    mass = 0.0
    for i in range(len(peptide) - 1):
        mass += AA_MASSES[peptide[i]]
        b_ions[f"b{i+1}"] = (mass + charge * MASS_PROTON) / charge
    mass = 0.0
    for i in range(len(peptide) - 1, 0, -1):
        mass += AA_MASSES[peptide[i]]
        y_ions[f"y{len(peptide)-i}"] = (mass + MASS_H2O + charge * MASS_PROTON) / charge
    return b_ions, y_ions

# Format spectrum data as text input for model
def spectrum_to_input_text(mz_list, intensity_list):
    pairs = sorted(zip(mz_list, intensity_list))
    lines = ["Experimental MS/MS spectrum (m/z, intensity):"]
    lines += [f"  m/z {mz:.4f}  intensity {intensity}" for mz, intensity in pairs]
    lines.append("Analyze this spectrum and identify the most likely peptide sequence.")
    return "\n".join(lines)

# Analyze spectrum and generate reasoning output
def generate_reasoning_output(peptide, mz_list, intensity_list, tolerance=0.05):
    b_ions, y_ions = theoretical_ions(peptide)
    matched_b, matched_y, isotopes, noise = [], [], [], []

    for mz, intensity in zip(mz_list, intensity_list):
        hit = False
        for name, tmz in b_ions.items():
            if abs(mz - tmz) <= tolerance:
                matched_b.append((name, mz, intensity)); hit = True; break
            elif abs(mz - (tmz + MASS_C13)) <= tolerance:
                isotopes.append((f"{name}+1", mz, intensity)); hit = True; break
        if hit:
            continue
        for name, tmz in y_ions.items():
            if abs(mz - tmz) <= tolerance:
                matched_y.append((name, mz, intensity)); hit = True; break
            elif abs(mz - (tmz + MASS_C13)) <= tolerance:
                isotopes.append((f"{name}+1", mz, intensity)); hit = True; break
        if not hit:
            noise.append((mz, intensity))

    top3 = " ".join(f"{mz:.4f}" for mz, _ in
                    sorted(zip(mz_list, intensity_list), key=lambda x: x[1], reverse=True)[:3])

    parts = [
        "Input Experimental MSMS Spectrum",
        f"Step 1 Spectrum summary. Number of peaks {len(mz_list)}.",
        f"mz range {min(mz_list):.2f} to {max(mz_list):.2f}." if mz_list else "mz range 0.00 to 0.00.",
        f"Most intense peaks {top3}.",
    ]

    if matched_b:
        parts += [f"{name} fragment consistent with substring {peptide[:int(name[1:])]} with a peak at {mz:.4f}."
                  for name, mz, _ in matched_b]
    else:
        parts.append("No b ion evidence found.")

    if matched_y:
        parts += [f"{name} fragment consistent with substring {peptide[-int(name[1:]):]} with a peak at {mz:.4f}."
                  for name, mz, _ in matched_y]
    else:
        parts.append("No y ion evidence found.")

    if isotopes:
        parts.append(f"Isotope pattern peaks {' '.join(f'{mz:.4f}' for _, mz, _ in isotopes)}.")
    else:
        parts.append("No clear isotopic satellite peaks detected.")

    if noise:
        sample = " ".join(f"{mz:.2f}" for mz, _ in sorted(noise)[:4])
        parts.append(f"Unassigned peaks include {sample} and {len(noise)} total.")
    else:
        parts.append("No unassigned peaks.")

    n_evidence = len(matched_b) + len(matched_y) + len(isotopes)
    confidence = "High" if n_evidence >= 5 else "Moderate" if n_evidence >= 3 else "Low"
    parts += [f"Evidence {n_evidence} supporting fragment peaks.", f"Confidence {confidence}.",
              f"Predicted peptide {peptide}."]

    paragraph = re.sub(r"\s+", " ", " ".join(parts)).strip()
    paragraph = re.sub(r"[^A-Za-z0-9\., ]+", " ", paragraph.replace("MS/MS", "MSMS").replace("m/z", "mz"))
    return paragraph.rstrip(".") + "."

# Check if peptide sequence is valid
def is_valid_peptide(peptide):
    return bool(peptide) and bool(re.fullmatch(r'[A-Z]+', peptide)) and all(aa in AA_MASSES for aa in peptide)

# Parse PARQUET file and extract spectrum data
def parse_parquet(parquet_path, max_spectra=100):

    pf = pq.ParquetFile(parquet_path)
    batch = next(pf.iter_batches(batch_size = max_spectra))
    df = batch.to_pandas()

    spectra = [
        {
            "id": r.PSMId,
            "peptide": r.peptide,
            "pepmass": r.precursor_mz,
            "charge": int(r.charge),
            "mz": r.mz_array.tolist(),
            "intensity": r.intensity_array.tolist(),
        }

        for r in df.itertuples()
    ]

    return spectra

# Process PARQUET file and save training dataset
def process_parquet(parquet_path, output_path, tolerance=0.05, verbose=True, num_examples=None):
    spectra = parse_parquet(parquet_path)
    dataset, skipped = [], 0

    for i, spec in enumerate(spectra):
        if num_examples and i >= num_examples:
            break
        peptide = spec.get("peptide")
        if not peptide or not is_valid_peptide(peptide):
            if verbose and peptide:
                print(f"  [SKIP] {spec['id']} – invalid peptide: '{peptide}'")
            skipped += 1
            continue
        dataset.append({
            "input":   spectrum_to_input_text(spec["mz"], spec["intensity"]),
            "output":  generate_reasoning_output(peptide, spec["mz"], spec["intensity"], tolerance=tolerance),
            "peptide": peptide,
            "id":   spec["id"],
        })
        if verbose:
            print(f"  [OK] {spec['id']} → {peptide}  ({len(spec['mz'])} peaks)")

    with open(output_path, "w") as f:
        for ex in dataset:
            f.write(json.dumps({"input": ex["input"], "output": ex["output"]}) + "\n")

    print(f"\nProcessed {len(dataset)}, skipped {skipped}. Saved → '{output_path}'")
    return dataset

if __name__ == "__main__":
    PARQUET_PATH    = "/disk/raptor-array/mshab008/msbenchmarking/hardness_index_0/filtered_spectra/hardness_index_0.parquet"
    OUTPUT_PATH = "spectra_training_data.jsonl"
    TOLERANCE   = 0.1
    VERBOSE     = True
    NUM_EXAMPLES = 10  # set to int to limit

    print(f"\nReading PARQUET: {PARQUET_PATH}")
    process_parquet(PARQUET_PATH, OUTPUT_PATH, tolerance=TOLERANCE, verbose=VERBOSE, num_examples=NUM_EXAMPLES)