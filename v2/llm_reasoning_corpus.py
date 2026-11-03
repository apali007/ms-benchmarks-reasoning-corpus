import os
import json
import re
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import argparse
from datetime import datetime

from constants import (
    AA_MASSES, NEUTRAL_LOSS_MASSES, MASS_PROTON,
    DATASET_METADATA_ROOT,
)

def extend_ladder(ion_type, matched_ion, label_peptide, precursor_mz, precursor_charge, explain_precursor):
    precursor_neutral_mass = precursor_mz * precursor_charge - precursor_charge * MASS_PROTON
    residues = list(reversed(label_peptide)) if ion_type == 'y' else list(label_peptide)
    n = len(residues)
    traces = []
    ions = []
    residues_found = []
    prev_mass = None

    # Position from N-terminus (1-based), so y and b refer to the same residue numbering
    def to_pos(i):
        return n - i + 1 if ion_type == 'y' else i

    for i in range(1, n + 1):
        peak_information = matched_ion.get(f"{ion_type}{i}")

        if peak_information is None:
            # The full length ion is never observed; recover the last residue from the precursor
            if i == n and prev_mass is not None:
                if explain_precursor:
                    precursor_text = (f"The neutral peptide mass is computed from the precursor m/z ({precursor_mz:.5f}) and charge "
                                      f"({precursor_charge}+): {precursor_mz:.5f} x {precursor_charge} - {precursor_charge} x "
                                      f"{MASS_PROTON:.5f} = {precursor_neutral_mass:.5f} Da. "
                                      )
                else:
                    precursor_text = ""
                if ion_type == 'y':
                    aa_mass = precursor_neutral_mass - prev_mass + MASS_PROTON
                    traces.append(
                        f"The remaining residue at the N-terminus comes from the precursor. {precursor_text}"
                        f"The neutral peptide mass ({precursor_neutral_mass:.5f} Da) minus y{n-1} ({prev_mass:.5f} Da) plus one proton "
                        f"({MASS_PROTON:.5f} Da) gives {aa_mass:.5f} Da, which matches the residue mass of {residues[i-1]}."
                    )
                else:
                    aa_mass = precursor_neutral_mass - NEUTRAL_LOSS_MASSES['H2O'] - prev_mass + MASS_PROTON
                    traces.append(
                        f"The remaining residue at the C-terminus comes from the precursor. {precursor_text}"
                        f"The neutral peptide mass ({precursor_neutral_mass:.5f} Da) minus H2O ({NEUTRAL_LOSS_MASSES['H2O']:.5f} Da), " 
                        f"minus b{n-1} ({prev_mass:.5f} Da), plus one proton ({MASS_PROTON:.5f} Da) gives {aa_mass:.5f} Da, "
                        f"which matches the residue mass of {residues[i-1]}."
                    )
                       
                residues_found.append({"pos": to_pos(i), "method": "precursor"})

            prev_mass = None
            continue

        (mz, intensity, charge, loss) = peak_information
        ions.append({"ion":f"{ion_type}{i}", "mz": round(float(mz), 5), "intensity": round(float(intensity), 5),  
                    "charge": charge, "loss": loss})
        
        ion_mass = mz * charge - (charge - 1) * MASS_PROTON
        if charge > 1:
            protons_text = f"{MASS_PROTON:.5f}" if charge == 2 else f"{charge-1} x {MASS_PROTON:.5f}"
            line = (f"{ion_type}{i} is observed at m/z {mz:.5f} as a {charge}+ ion. "
                    f"Converting it to a singly charged ion gives {mz:.5f} x {charge} - {protons_text} = {ion_mass:.5f} Da.")
        else:
            line = f"{ion_type}{i} is observed at m/z {mz:.5f}."
        
        if loss is not None:
            loss_mass = NEUTRAL_LOSS_MASSES[loss]
            line += (f" This peak has lost {loss} ({loss_mass:.5f} Da), "
                     f"so adding it back gives {ion_mass + loss_mass:.5f} Da.")
            ion_mass +=loss_mass
            
        if i == 1: 
            if ion_type == 'y':
                aa_mass = ion_mass - NEUTRAL_LOSS_MASSES['H2O'] - MASS_PROTON
                line += (f" Since y1 contains only the last residue plus H2O and a proton, subtracting H2O "
                f"({NEUTRAL_LOSS_MASSES['H2O']:.5f} Da) and a proton ({MASS_PROTON:.5f} Da) leaves {aa_mass:.5f} Da, "
                f"which matches the residue mass of {residues[0]}.")
            else:
                aa_mass = ion_mass - MASS_PROTON
                line += (f" Since b1 contains only the first residue plus a proton, subtracting the proton "
                f"({MASS_PROTON:.5f} Da) leaves {aa_mass:.5f} Da, which matches the residue mass of {residues[0]}.")
            residues_found.append({"pos": to_pos(1), "method": "first"})

        elif prev_mass is not None:
            line += (f" The mass difference between {ion_type}{i-1} and {ion_type}{i} is {ion_mass:.5f} - {prev_mass:.5f} = "
            f"{ion_mass - prev_mass:.5f} Da, which matches the residue mass of {residues[i-1]}.")
            residues_found.append({"pos": to_pos(i), "method": "diff"})

        traces.append(line)
        prev_mass = ion_mass

    return traces, ions, residues_found

def derive_peptide(matched_y, matched_b, label_peptide, precursor_mz, precursor_charge):
    y_traces, y_ions, y_residues = extend_ladder('y', matched_y, label_peptide, precursor_mz, precursor_charge, explain_precursor=True)

    # Explain the precursor conversion in the b ladder only if the y ladder did not already use it.   
    y_used_precursor = any(r['method'] == 'precursor' for r in y_residues)
    b_traces, b_ions, b_residues = extend_ladder('b', matched_b, label_peptide, precursor_mz, precursor_charge, explain_precursor= not y_used_precursor)

    intro =(f"To determine the peptide sequence, I derive residues from the y-ion ladder, starting at the C-terminus, "
            "and then the b-ion ladder, starting at the N-terminus. Each residue is identified from the mass " 
            f"difference between consecutive ions.")

    if not y_traces:
        y_traces = ["No y ions were found in this spectrum."]
    if not b_traces:
        b_traces = ["No b ions were found in this spectrum."]

    final_sequence = ''.join(label_peptide)
    reasoning_trace = (
        [intro, "", "Deriving residues from the y-ion ladder (C-terminus to N-terminus):"] + y_traces +
        ["", "Deriving residues from the b-ion ladder (N-terminus to C-terminus):"] + b_traces +
        ["", f"Combining the y-ion and b-ion ladders gives the final sequence {final_sequence}."]
    )

    ladder_info = {
        "y_ions": y_ions,
        "b_ions": b_ions,
        "y_residues": y_residues,
        "b_residues": b_residues,
    }

    return label_peptide, reasoning_trace, ladder_info
    
def spectrum_to_input_text(mz_list, intensity_list, precursor_mz, precursor_charge):
    unique_peaks = sorted(set(zip(mz_list, intensity_list)), key=lambda p: p[0])
    # peaks = sorted(zip(mz_list, intensity_list), key=lambda p: p[0])
    lines = [f"Precursor mz: {precursor_mz:.5f} Da, Precursor Charge: {precursor_charge}"]
    lines += ["Experimental MS/MS spectrum (m/z, intensity):"]
    lines += [f"m/z {mz:.5f}, intensity {intensity:.5f}" for mz, intensity in unique_peaks]
    lines.append("Analyze this spectrum and identify the most likely peptide sequence.")
    return "\n".join(lines)

def is_valid_peptide(peptide):
    return bool(peptide) and bool(re.fullmatch(r'[A-Z]+', peptide)) and all(aa in AA_MASSES for aa in peptide)

# Compute theoretical b and y ions for a peptide sequence
def theoretical_ions(peptide, aa_masses, charge=1):
    b_ions, y_ions = {}, {}
    mass = 0.0
    for i in range(len(peptide) - 1):
        mass += aa_masses[peptide[i]]
        b_ions[f"b{i+1}"] = (mass + charge * MASS_PROTON) / charge
    mass = 0.0
    for i in range(len(peptide) - 1, 0, -1):
        mass += aa_masses[peptide[i]]
        y_ions[f"y{len(peptide)-i}"] = (mass + NEUTRAL_LOSS_MASSES['H2O'] + charge * MASS_PROTON) / charge
    return b_ions, y_ions

def load_search_params(dataset):
    metadata_path = os.path.join(DATASET_METADATA_ROOT, f"{dataset}.json")
    with open(metadata_path) as f:
        metadata = json.load(f)

    aa_masses = dict(AA_MASSES)
    for entry in metadata.get('modifications', '').split(','):
        entry = entry.strip()
        if not entry:
            continue
        m = re.match(r'^([A-Z]+)([+-][\d.]+)$', entry)
        if m is None:
            raise ValueError(f"Cannot parse modification '{entry}' in {metadata_path}")
        residues, delta = m.groups()
        for aa in residues:
            aa_masses[aa] += float(delta)


    return {
        'aa_masses': aa_masses,
        'fragment_tolerance': float(metadata['Fragment_Tolerance']),
        'precursor_window': float(metadata['Precursor_mz_Tolerance']),
        'precursor_window_type': metadata['Precursor_mz_Unit'].strip().lower(),
    }

# Parse PARQUET file and extract spectrum data
def parse_parquet(parquet_path, num_examples=100):

    columns = ['dataset', 'PSMId', 'peptide', 'modified_peptide', 'precursor_mz', 'charge', 'mz_array', 'intensity_array']
    pf = pq.ParquetFile(parquet_path)
    spectra = []
    n_selected = 0
    n_modified = 0

    for batch in pf.iter_batches(batch_size=10000, columns=columns):

        df = batch.to_pandas()

        for r in df.itertuples():
            if r.peptide != r.modified_peptide:
                n_modified += 1
                continue

            spectra.append({
                "dataset": r.dataset,
                "psmid": r.PSMId,
                "peptide": r.peptide,
                "precursor_mz": r.precursor_mz,
                "precursor_charge": int(r.charge),
                "mz": r.mz_array.tolist(),
                "intensity": r.intensity_array.tolist(),
            })
            
            n_selected += 1

            if num_examples is not None and n_selected == num_examples:
                print(f"Selected {n_selected}, skipped {n_modified} modified peptides")
                return spectra

    print(f"Selected {n_selected}, skipped {n_modified} modified peptides")
    return spectra

def match_ions(peptide, mz_list, intensity_list, aa_masses, fragment_tolerance, max_fragment_charge):
    # Each theoretical ion: (loss_rank, name, charge, loss, m/z). Intact ions have loss_rank 0, neutral losses 1.
    theoretical = []
    peaks = list(zip(mz_list, intensity_list))

    for charge in range(1, max_fragment_charge + 1):
        b_ions, y_ions = theoretical_ions(peptide, aa_masses, charge)
        for name, tmz in {**b_ions, **y_ions}.items():
            theoretical.append((0, name, charge, None, tmz))
            for loss, loss_mass in NEUTRAL_LOSS_MASSES.items():
                theoretical.append((1, name, charge, loss, tmz - loss_mass / charge))

    # All (ion, charge, peak) pairs within tolerance
    candidates = []
    for loss_rank, name, charge, loss, tmz in theoretical:
        for peak_idx, (mz, intensity) in enumerate(peaks):
            error = abs(mz-tmz)
            if error <= fragment_tolerance:
                candidates.append((loss_rank, charge, error, -intensity, name, loss, peak_idx))

    # Intact before neutral loss, lower charge first, then closest, then most intense; 
    # each ion and each peak used at most once
    candidates.sort()
    assigned, assigned_peaks = {}, set()
    for loss_rank, charge, error, _, name, loss, peak_idx in candidates:
        if (name, charge, loss) in assigned or peak_idx in assigned_peaks:
            continue
        assigned[(name, charge, loss)] = peaks[peak_idx]
        assigned_peaks.add(peak_idx)

    # One peak per ladder position: prefer intact over neutral loss, then the lowest charge
    matched_b, matched_y = {}, {}
    for (name, charge, loss), (mz, intensity) in sorted(assigned.items(), key=lambda kv: (kv[0][2] is not None, kv[0][1])):
        target = matched_b if name.startswith('b') else matched_y
        if name not in target:
            target[name] = (mz, intensity, charge, loss)

    noise = [peak for idx, peak in enumerate(peaks) if idx not in assigned_peaks]

    return matched_y, matched_b, assigned, noise

def normalize(s):
    return s.replace('I', 'L')

def process_input_data(input_data_path, output_path, num_examples, save_dataset):
    spectra = parse_parquet(input_data_path, num_examples)
    records = []
    params_by_dataset = {}

    for total, spec in enumerate(spectra, start=1):
        peptide = spec['peptide']
        dataset_name = spec['dataset']
        if dataset_name not in params_by_dataset:
            params_by_dataset[dataset_name] = load_search_params(dataset_name)
        params = params_by_dataset[dataset_name]

        matched_y, matched_b, annotations, noise = match_ions(
            peptide, spec['mz'], spec['intensity'], 
            params['aa_masses'], params['fragment_tolerance'], spec['precursor_charge'])
        
        final_sequence, reasoning_trace, ladder_info = derive_peptide(
            matched_y,
            matched_b,
            peptide,
            spec['precursor_mz'],
            spec['precursor_charge'],
        )

        if VERBOSE:
            print(f"\n==== {spec['psmid']} | {peptide} ====")
            print('\n'.join(reasoning_trace))

        if save_dataset:
            records.append({
            'id': spec['psmid'],
            'dataset': dataset_name,
            'input': spectrum_to_input_text(
                spec['mz'], 
                spec['intensity'], 
                spec['precursor_mz'], 
                spec['precursor_charge'],),
            'output': '\n'.join(reasoning_trace),
            'peptide': peptide,
            **ladder_info,
            })

            if total % 1000 == 0:
                with open(output_path, 'w') as f:
                    json.dump(records, f, indent=2)
                print(f"Progress: {total} processed")

    print(f"\nProcessed {len(spectra)} spectra")
    if save_dataset:
        with open(output_path, 'w') as f:
            json.dump(records, f, indent=2)
        print(f"Saved {len(records)} training examples to {output_path}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate the reasoning corpus for LLM to perform de novo peptide sequencing")
    parser.add_argument("--input_data_path",
                         default="/disk/raptor-array/mshab008/msbenchmarking/hardness_index_0/filtered_spectra/hardness_index_0.parquet")
    parser.add_argument("--output_dir",
                         default="v2/output")
    parser.add_argument("--num_examples", type=int, default=10,
                         help="Number of spectra to keep in the output dataset.")
    parser.add_argument("--no_save", action="store_true",
                         help="Process spectra without writing a dataset file (stats-only run).")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    input_data_path = args.input_data_path
    VERBOSE = args.verbose
    NUM_EXAMPLES = args.num_examples if args.num_examples > 0 else None
    SAVE_DATASET = not args.no_save

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    n_tag = "all" if NUM_EXAMPLES is None else NUM_EXAMPLES
    os.makedirs(args.output_dir, exist_ok=True)
    OUTPUT_PATH = os.path.join(args.output_dir, f"{timestamp}_{n_tag}_llm_reasoning_corpus.json")

    process_input_data(input_data_path, OUTPUT_PATH, NUM_EXAMPLES, SAVE_DATASET)