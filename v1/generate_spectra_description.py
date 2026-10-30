import json
import re
import numpy as np


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

# Generate synthetic mass spectrum with peaks and noise
def build_synthetic_spectrum(peptide, b_coverage=0.6, y_coverage=0.7,
                             isotope_prob=0.4, noise_peaks=4,
                             intensity_noise_frac=0.05, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    b_ions, y_ions = theoretical_ions(peptide)
    mz_list, int_list = [], []

    def add(mz, base_intensity):
        jitter = rng.normal(0, intensity_noise_frac * base_intensity)
        mz_list.append(round(mz, 4))
        int_list.append(max(50, int(base_intensity + jitter)))

    for mz in b_ions.values():
        if rng.random() < b_coverage:
            intensity = rng.integers(600, 2000)
            add(mz, intensity)
            if rng.random() < isotope_prob:
                add(mz + MASS_C13, intensity * rng.uniform(0.1, 0.35))

    for mz in y_ions.values():
        if rng.random() < y_coverage:
            intensity = rng.integers(800, 2500)
            add(mz, intensity)
            if rng.random() < isotope_prob:
                add(mz + MASS_C13, intensity * rng.uniform(0.1, 0.35))

    all_mz = mz_list or [100.0]
    for _ in range(noise_peaks):
        add(round(rng.uniform(min(all_mz) * 0.8, max(all_mz) * 1.1), 2),
            rng.integers(40, 200))

    return mz_list, int_list

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

# Create training dataset from peptides
def generate_dataset(peptides, samples_per_peptide=3, seed=42):
    rng = np.random.default_rng(seed)
    dataset = []
    for peptide in peptides:
        for _ in range(samples_per_peptide):
            mz_list, int_list = build_synthetic_spectrum(
                peptide,
                b_coverage=rng.uniform(0.4, 0.85),
                y_coverage=rng.uniform(0.5, 0.9),
                isotope_prob=rng.uniform(0.2, 0.6),
                noise_peaks=rng.integers(2, 7),
                rng=rng,
            )
            dataset.append({
                "input":   spectrum_to_input_text(mz_list, int_list),
                "output":  generate_reasoning_output(peptide, mz_list, int_list),
                "peptide": peptide,
            })
    return dataset

# Check if peptide sequence is valid
def is_valid_peptide(peptide):
    return bool(peptide) and bool(re.fullmatch(r'[A-Z]+', peptide)) and all(aa in AA_MASSES for aa in peptide)

# Parse MGF file and extract spectrum data
def parse_mgf(mgf_path, max_spectra=100):
    spectra, current, count = [], None, 0
    with open(mgf_path) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line == "BEGIN IONS":
                count += 1
                if count > max_spectra:
                    break
                current = {"title": "", "peptide": None, "charge": 1, "pepmass": 0.0, "mz": [], "intensity": []}
            elif line == "END IONS":
                if current is not None:
                    spectra.append(current)
                current = None
            elif current is None:
                continue
            elif line.startswith("TITLE="):
                current["title"] = line[6:]
            elif line.startswith("SEQ="):
                current["peptide"] = ''.join(filter(str.isupper, line[4:]))
            elif line.startswith("PEPMASS="):
                current["pepmass"] = float(line.split("=", 1)[1].split()[0])
            elif line.startswith("CHARGE="):
                try:
                    current["charge"] = int(line.split("=", 1)[1].strip().rstrip("+-"))
                except ValueError:
                    current["charge"] = 1
            elif line[0].isdigit() or line[0] == ".":
                parts = line.split()
                if len(parts) >= 2:
                    current["mz"].append(float(parts[0]))
                    current["intensity"].append(float(parts[1]))
    return spectra

# Process MGF file and save training dataset
def process_mgf(mgf_path, output_path, tolerance=0.05, verbose=True, num_examples=None):
    spectra = parse_mgf(mgf_path)
    dataset, skipped = [], 0

    for i, spec in enumerate(spectra):
        if num_examples and i >= num_examples:
            break
        peptide = spec.get("peptide")
        if not peptide or not is_valid_peptide(peptide):
            if verbose and peptide:
                print(f"  [SKIP] {spec['title']} – invalid peptide: '{peptide}'")
            skipped += 1
            continue
        dataset.append({
            "input":   spectrum_to_input_text(spec["mz"], spec["intensity"]),
            "output":  generate_reasoning_output(peptide, spec["mz"], spec["intensity"], tolerance=tolerance),
            "peptide": peptide,
            "title":   spec["title"],
        })
        if verbose:
            print(f"  [OK] {spec['title']} → {peptide}  ({len(spec['mz'])} peaks)")

    with open(output_path, "w") as f:
        for ex in dataset:
            f.write(json.dumps({"input": ex["input"], "output": ex["output"]}) + "\n")

    print(f"\nProcessed {len(dataset)}, skipped {skipped}. Saved → '{output_path}'")
    return dataset


if __name__ == "__main__":
    MGF_PATH    = "/disk/raptor-array/SaeedLab-Data/Mass-Spec/SpectralLibraries/Massive-KB-v1/LIBRARY_CREATION_AUGMENT_LIBRARY_TEST-82c0124b-download_filtered_mgf_library-main.mgf"
    OUTPUT_PATH = "spectra_training_data.jsonl"
    TOLERANCE   = 0.1
    VERBOSE     = False
    NUM_EXAMPLES = None  # set to int to limit

    print(f"\nReading MGF: {MGF_PATH}")
    process_mgf(MGF_PATH, OUTPUT_PATH, tolerance=TOLERANCE, verbose=VERBOSE, num_examples=NUM_EXAMPLES)