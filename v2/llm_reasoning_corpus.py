import os
import json
import re
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import argparse

from constants import (
    AA_MASSES, NEUTRAL_LOSS_MASSES, MASS_PROTON,
    DATASET_METADATA_ROOT,
    CASE_CLEAN_B_ION_LADDER, CASE_CLEAN_Y_ION_LADDER,
    CASE_COMPLEMENTARY_GAPS, CASE_UNRESOLVABLE_MULTI_GAP,
    CASE_SINGLE_GAP_RESOLVED, CASE_UNRESOLVABLE_SINGLE_GAP_MASS_MISMATCH,
    CASES
)


def ordinal(n):
    #Convert an integer to its ordinal string: 1 -> '1st', 2 -> '2nd', 11 -> '11th', etc.
    if 11 <= (n % 100) <= 13:
        suffix = 'th'
    else:
        suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')
    return f"{n}{suffix}"

def format_ladder_sequence_dash(derived_aa):
    seq_str = '-'.join(aa if aa is not None else '?' for aa in derived_aa)
    return seq_str

def format_ladder_sequence(derived_aa):
    seq_str = ''.join(aa if aa is not None else '?' for aa in derived_aa)
    return seq_str

def handle_terminal_ion(ion_type, n, prev_mz, precursor_neutral_mass, predicted_aa_mass, terminal_aa):
    
    terminal_ion_trace = (f"No fragment ion was found for {ion_type}{n}. Since this is the terminal residue, "
        f"it is instead recovered from {ion_type}{n-1} and the precursor mass: precursor-derived neutral mass {precursor_neutral_mass:.5f} Da, "
        f"minus the {ion_type}{n-1} fragment mass ({prev_mz:.5f} Da), plus one proton ({MASS_PROTON:.5f} Da), gives {predicted_aa_mass:.5f} Da. This identifies the residue mass as {terminal_aa}."
        )
    if VERBOSE:
        print(terminal_ion_trace)
    
    return terminal_ion_trace

def validate_precursor(final_sequence, precursor_neutral_mass, aa_masses, precursor_tol_da):

    theoretical_mass_with_water = sum(aa_masses[r] for r in final_sequence) +  NEUTRAL_LOSS_MASSES['H2O']
    diff = abs(precursor_neutral_mass - theoretical_mass_with_water)
    is_valid = diff <= precursor_tol_da

    return is_valid, diff, theoretical_mass_with_water

def final_precursor_check(final_sequence, precursor_neutral_mass, aa_masses, precursor_tol_da):

    is_precursor_valid, diff, theoretical_mass_with_water = validate_precursor(final_sequence, precursor_neutral_mass, aa_masses, precursor_tol_da)

    if not is_precursor_valid:
        return "", False
    else:
        trace = (f"Validating against the precursor mass: Theoretical mass of {final_sequence} (sum of residue masses plus one water molecule): {theoretical_mass_with_water:.5f} Da. "
        f"Precursor-derived neutral mass: {precursor_neutral_mass:.5f} Da. Difference: {diff:.5f} Da is within the {precursor_tol_da:.5f} Da tolerance, confirming the derived sequence is consistent with the measured precursor mass. ")

        return trace, True
    
def handle_ion_gaps(b_derived, y_derived, label_peptide, n):
    trace = []
    for i in range(n):
        if b_derived[i] is None:
            if y_derived[i] is not None:
                b_derived[i] = y_derived[i]
                trace.append(f"The residue at position {ordinal(i+1)} (from the N-terminus) was identified by the y-ion ladder as {b_derived[i]}. ")
                if VERBOSE:
                    print(f"The residue at position {ordinal(i+1)} (from the N-terminus) was identified by the y-ion ladder as {b_derived[i]}. ")
            else:
                trace.append(f"The residue at position {ordinal(i+1)} (from the N-terminus) cannot be resolved by the y-ion ladder. ")
                if VERBOSE:
                    print(f"The residue at position {ordinal(i+1)} (from the N-terminus) cannot be resolved by the y-ion ladder. ")

    b_seq_str = format_ladder_sequence_dash(b_derived)
    ion_gap_trace = '\n'.join(trace)

    return b_seq_str, b_derived, ion_gap_trace

def handle_first_position(mz, ion_type, ladder_residues):
    if ion_type == 'y':
        predicted_aa_mass = mz-NEUTRAL_LOSS_MASSES['H2O']-MASS_PROTON
        trace = (
                f"Subtracting H2O ({NEUTRAL_LOSS_MASSES['H2O']:.5f} Da) and a proton ({MASS_PROTON:.5f} Da) "
                f"gives {predicted_aa_mass:.5f} Da. "
                f"This identifies the 1st residue from the C-terminus as {ladder_residues[0]}."
            )
    else:
        predicted_aa_mass = mz-MASS_PROTON
        trace = (
        f"Subtracting a proton ({MASS_PROTON:.5f} Da) gives {mz-MASS_PROTON:.5f} Da. "
        f"This identifies the 1st residue from the N-terminus as {ladder_residues[0]}."
        )
    
    return trace

def extend_ladder(ion_type, matched_ion, label_peptide, precursor_neutral_mass, precursor_charge, aa_masses, precursor_tol_da):
    residues = list(label_peptide)
    n = len(residues)
    if ion_type == 'y':
        ladder_residues = list(reversed(residues))  
        terminal = 'C'
    else: 
        ladder_residues = residues
        terminal = 'N'

    traces = []

    # Initializing the position of fragment ions to be derived with none.
    derived_aa_pos = [None] * n
    prev_mz = None

    for i in range(1, n + 1):
        key = f"{ion_type}{i}"
        peak_information = matched_ion.get(key)
        if VERBOSE:
                print(f"{ion_type}{i}, peak_info: {peak_information}")
        
        if peak_information is None:
            if i == n and prev_mz is not None:
                if ion_type == 'y':
                    predicted_aa_mass = precursor_neutral_mass - prev_mz + MASS_PROTON
                else:
                    predicted_aa_mass = precursor_neutral_mass - NEUTRAL_LOSS_MASSES['H2O'] - prev_mz + MASS_PROTON
                is_precursor_validated, _ , _ = validate_precursor(label_peptide, precursor_neutral_mass, aa_masses, precursor_tol_da)                
                if is_precursor_validated:
                    traces.append(handle_terminal_ion(ion_type, n, prev_mz, precursor_neutral_mass, predicted_aa_mass, ladder_residues[i-1]))
                    derived_aa_pos[i-1] = ladder_residues[i-1] 
                else:
                    traces.append(f"No fragment ion was found for {ion_type}{i}. ")
            else:
                traces.append(f"No fragment ion was found for {ion_type}{i}. ")

            prev_mz = None
            continue

        (mz, intensity) = peak_information
        if i == 1: 
            traces.append(handle_first_position(mz, ion_type, ladder_residues))
            derived_aa_pos[i-1] = ladder_residues[i-1]
        else:
            if prev_mz is not None:
                predicted_aa_mass = mz-prev_mz
                traces.append(
                        f"The difference between {ion_type}{i} and {ion_type}{i-1} ({mz:.5f} - {prev_mz:.5f} = {predicted_aa_mass:.5f} Da) identifies the {ordinal(i)} residue from the {terminal}-terminus as {ladder_residues[i-1]}. "
                    )
                derived_aa_pos[i-1] = ladder_residues[i-1]
            # else:
            #     # assign amino acid if prev_mz is none but the peak exists in the labeled dataset. But you cannot compute the difference ss, you just use the labeled peptide position.
            #     traces.append(
            #     f"This identifies the {ordinal(i)} residue from the {terminal}-terminus as {ladder_residues[i-1]}. "
            #         )
            #     derived_aa_pos[i-1] = ladder_residues[i-1]


        prev_mz = mz

    if VERBOSE:
        detected_positions = [f"{ion_type}{i+1}" for i, aa in enumerate(derived_aa_pos) if aa is not None]
        gap_positions = [f"{ion_type}{i+1}" for i, aa in enumerate(derived_aa_pos) if aa is None]
        print(f"detected_positions: {', '.join(detected_positions)}")
        print(f"gap_positions: {', '.join(gap_positions)}")

    derived_aa_pos = list(reversed(derived_aa_pos)) if ion_type == 'y' else derived_aa_pos  
        
    return derived_aa_pos, traces

def derive_peptide(matched_y, matched_b, label_peptide, precursor_mz, precursor_charge, params):
    label_residues = list(label_peptide)
    precursor_neutral_mass = precursor_mz * precursor_charge - precursor_charge * MASS_PROTON
    aa_masses = params['aa_masses']
    fragment_tolerance = params['fragment_tolerance']
    precursor_tol_da =  precursor_tolerance_da(params, precursor_neutral_mass, precursor_charge)
    n = len(label_peptide)

    if VERBOSE:
        print(f"Length of the peptide: {n}")
        print(f"precursor_mz: {precursor_mz:.5f}")
        print(f"precursor_charge: {precursor_charge}")
        print(f"precursor_tolerance: {precursor_tol_da:.5f} Da")
        print(f"fragment_tolerance: {fragment_tolerance:.5f} Da")
        
    reasoning_trace = []

    tolerance_trace = (f"Analyzing spectra for the y-ion and b-ion ladder to predict the peptide sequence, "
        f"using a fragment mass tolerance of {fragment_tolerance:.5f} Da and a precursor mass tolerance of "
        f"{params['precursor_window']:g} {params['precursor_window_type']} ({precursor_tol_da:.5f} Da).")
    reasoning_trace.append(tolerance_trace)

    if VERBOSE:
        print(tolerance_trace)

    reasoning_trace.append("First, trying to derive the sequence using the y-ion ladder (C-terminus to N-terminus): ")
    y_derived, y_traces = extend_ladder('y', matched_y, label_peptide, precursor_neutral_mass, precursor_charge, aa_masses, precursor_tol_da)
    reasoning_trace.append('\n'.join(y_traces))

    y_seq_str = format_ladder_sequence_dash(y_derived)
    y_has_gap = any(aa is None for aa in y_derived)

    reasoning_trace.append(f"The y-ion ladder identifies this sequence (N-terminus to C-terminus): {y_seq_str}. ")

    if not y_has_gap:
        final_sequence = "".join(y_derived)
        precursor_validation_trace, is_precursor_match = final_precursor_check(final_sequence, precursor_neutral_mass, aa_masses, precursor_tol_da)
        reasoning_trace.append(precursor_validation_trace)
        if is_precursor_match:
            reasoning_trace.append(f"No need to derive a b-ion ladder as the complete sequence is confirmed. Final predicted sequence: {final_sequence}. ")
        is_match = (label_peptide == final_sequence)
        return is_match, is_precursor_match, final_sequence, reasoning_trace, CASE_CLEAN_Y_ION_LADDER

    reasoning_trace.append("Next, trying to derive the sequence using the b-ion ladder (N-terminus to C-terminus): ")
    b_derived, b_traces = extend_ladder('b', matched_b, label_peptide, precursor_neutral_mass, precursor_charge, aa_masses, precursor_tol_da)
    reasoning_trace.append('\n'.join(b_traces))

    b_seq_str = format_ladder_sequence_dash(b_derived)
    b_has_gap = any(aa is None for aa in b_derived)

    reasoning_trace.append(f"The b-ion ladder identifies this sequence (N-terminus to C-terminus): {b_seq_str}. ")
    
    if not b_has_gap:
        final_sequence = "".join(b_derived)
        precursor_validation_trace, is_precursor_match = final_precursor_check(final_sequence, precursor_neutral_mass, aa_masses, precursor_tol_da)
        reasoning_trace.append(precursor_validation_trace)
        if is_precursor_match:
            reasoning_trace.append(f"Final predicted sequence: {final_sequence}. ")
        is_match = (label_peptide == final_sequence)
        return is_match, is_precursor_match, final_sequence, reasoning_trace, CASE_CLEAN_B_ION_LADDER

    reasoning_trace.append("Resolving the remaining b-ion ladder gaps using the y-ion ladder:")
    b_seq_str_after_handle, b_derived_after_handle, ion_gap_trace = handle_ion_gaps(b_derived, y_derived, label_peptide, n)
    reasoning_trace.append(ion_gap_trace)

    gap_indices = [i for i, aa in enumerate(b_derived_after_handle) if aa is None]
    gap_count = len(gap_indices)    

    if gap_count == 0:
        final_sequence_after_handle = "".join(b_derived_after_handle)
        reasoning_trace.append(f"Combining both ladders gives the full predicted sequence: {final_sequence_after_handle}. ")
        precursor_validation_trace, is_precursor_match = final_precursor_check(final_sequence_after_handle, precursor_neutral_mass, aa_masses, precursor_tol_da)
        reasoning_trace.append(precursor_validation_trace)
        if is_precursor_match:
            reasoning_trace.append(f"Final predicted sequence: {final_sequence_after_handle}. ")
        is_match = (label_peptide == final_sequence_after_handle)
        return is_match, is_precursor_match, final_sequence_after_handle, reasoning_trace, CASE_COMPLEMENTARY_GAPS
    
    if gap_count == 1:
        gap_pos = gap_indices[0]
        aa_gap = label_residues[gap_pos]
        sum_identified_residues = sum(aa_masses[aa] for aa in b_derived_after_handle if aa is not None)
        identified_residues = '+'.join(aa for aa in b_derived_after_handle if aa is not None)
        predicted_aa_mass = precursor_neutral_mass - sum_identified_residues - NEUTRAL_LOSS_MASSES['H2O']
        is_precursor_validated, _ , _ = validate_precursor(label_peptide, precursor_neutral_mass, aa_masses, precursor_tol_da)
        if is_precursor_validated:
            incomplete_sequence = format_ladder_sequence(b_derived_after_handle) 
            reasoning_trace.append(
            f"Combining both ladders gives the incomplete predicted sequence: {incomplete_sequence}. "
            f"Since exactly one residue remains unresolved, it can be recovered directly from the precursor mass: "
            f"precursor-derived neutral mass {precursor_neutral_mass:.5f} Da, minus water ({NEUTRAL_LOSS_MASSES['H2O']:.5f} Da), minus the sum of the already-identified residues ({identified_residues} = {sum_identified_residues:.5f} Da), "
            f"leaves {predicted_aa_mass:.5f} Da, matching the mass of {aa_gap}. This identifies the residue at position {gap_pos+1} as {aa_gap}."
            )

            b_derived_after_handle[gap_pos] = aa_gap 
            final_sequence_after_handle = "".join(b_derived_after_handle)

            reasoning_trace.append(f"Combining all resolved positions gives the complete predicted sequence: {final_sequence_after_handle}. ")

            precursor_validation_trace, is_precursor_match = final_precursor_check(final_sequence_after_handle, precursor_neutral_mass, aa_masses, precursor_tol_da)
            reasoning_trace.append(precursor_validation_trace)
            if is_precursor_match:
                reasoning_trace.append(f"Final predicted sequence: {final_sequence_after_handle}. ")
            is_match = (label_peptide == final_sequence_after_handle)
            return is_match, is_precursor_match, final_sequence_after_handle, reasoning_trace, CASE_SINGLE_GAP_RESOLVED
            
        else:
            final_sequence_after_handle = format_ladder_sequence(b_derived_after_handle) 
            reasoning_trace.append(f"Combining both ladders gives the incomplete predicted sequence: {final_sequence_after_handle}. The complete sequence cannot be predicted. ")
            reasoning_trace.append(f"No final sequence is reported for this spectrum. ")
            is_match = (label_peptide == final_sequence_after_handle)
            return is_match, False, final_sequence_after_handle, reasoning_trace, CASE_UNRESOLVABLE_SINGLE_GAP_MASS_MISMATCH

    else:
        
        final_sequence_after_handle = format_ladder_sequence(b_derived_after_handle) 
        reasoning_trace.append(f"Combining both ladders gives the incomplete predicted sequence: {final_sequence_after_handle}. The complete sequence cannot be predicted. ")
        reasoning_trace.append(f"No final sequence is reported for this spectrum. ")
        is_match = (label_peptide == final_sequence_after_handle)
        return is_match, False, final_sequence_after_handle, reasoning_trace, CASE_UNRESOLVABLE_MULTI_GAP
    
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


def precursor_tolerance_da(params, precursor_neutral_mass, precursor_charge):
    window = params['precursor_window']
    window_type = params['precursor_window_type']
    if window_type == 'ppm':
        return window * precursor_neutral_mass / 1e6
    if window_type == 'da':
        return window
    raise ValueError(f"Unknown precursor tolerance unit: {window_type}")



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

def match_ions(peptide, mz_list, intensity_list, aa_masses, fragment_tolerance):
    b_ions, y_ions = theoretical_ions(peptide, aa_masses)
    theoretical = {**b_ions, **y_ions}
    peaks = list(zip(mz_list, intensity_list))

    # All (ion, peak) pairs within tolerance
    candidates = []
    for name, tmz in theoretical.items():
        for peak_idx, (mz, intensity) in enumerate(peaks):
            error = abs(mz-tmz)
            if error <= fragment_tolerance:
                candidates.append((error, -intensity, name, peak_idx))

    # Assign the closest pairs first; each ion and each peak used at most once
    candidates.sort()
    assigned_ions, assigned_peaks = {}, set()
    for error, _, name, peak_idx in candidates:
        if name in assigned_ions or peak_idx in assigned_peaks:
            continue
        assigned_ions[name] = peaks[peak_idx]
        assigned_peaks.add(peak_idx)

    matched_b = {name: peak for name, peak in assigned_ions.items() if name.startswith('b')}
    matched_y = {name: peak for name, peak in assigned_ions.items() if name.startswith('y')}
    noise = [peak for idx, peak in enumerate(peaks) if idx not in assigned_peaks]

    return matched_y, matched_b, noise

def normalize(s):
    return s.replace('I', 'L')

def process_input_data(input_data_path, output_path, num_examples, save_dataset):
    spectra = parse_parquet(input_data_path, num_examples)

    total = 0
    count_label_match = 0
    count_precursor_match = 0
    dataset = []
    case_counts = {}
    params_by_dataset = {}
    n_not_saved = 0

    for spec in spectra:
        total += 1
        peptide = spec['peptide']
        mz_list = spec['mz']
        intensity_list = spec['intensity']
        dataset_name = spec['dataset']
        if dataset_name not in params_by_dataset:
            params_by_dataset[dataset_name] = load_search_params(dataset_name)
        params = params_by_dataset[dataset_name]

        if VERBOSE:
            print(f"\n=====NEW PEPTIDE=====")
            print(f"label Peptide: {peptide}")
            print(f"PSM id: {spec['psmid']}")
    
        matched_y, matched_b, noise = match_ions(peptide, mz_list, intensity_list, params['aa_masses'], params['fragment_tolerance'])
        
        is_label_match, is_precursor_match, predicted_sequence, reasoning_trace, case_label = derive_peptide(
            matched_y,
            matched_b,
            spec['peptide'],
            spec['precursor_mz'],
            spec['precursor_charge'],
            params,
        )

        complete_resoning_trace = '\n'.join(reasoning_trace)

        if VERBOSE:
            print(f"is_label_match: {is_label_match}")
            print(f"is_precursor_match: {is_precursor_match}")
            print(f"predicted_sequence: {predicted_sequence}")
            print(f"case_label: {case_label}")
        

        if is_label_match: 
            count_label_match += 1
            if VERBOSE:
                print(f"predicted sequence is correctly matches with the label peptide")

        if is_precursor_match: 
            count_precursor_match += 1
            if VERBOSE:
                print(f"predicted sequence is correctly matches with the precursor information")

        case_counts[case_label] = case_counts.get(case_label, 0) + 1

        if save_dataset:
 
            dataset.append({
            'id': spec['psmid'],
            'input': spectrum_to_input_text(
                spec['mz'], 
                spec['intensity'], 
                spec['precursor_mz'], 
                spec['precursor_charge'],),
            'output': complete_resoning_trace,
            'peptide': peptide,
            'pred_seq': predicted_sequence,
            'is_label_match': is_label_match,
            'is_precursor_match': is_precursor_match,
            'case_label': case_label,
            })
        
            # Store correct peptides
            if total % 1000 == 0:
                with open(output_path, 'w') as f:
                    json.dump(dataset, f, indent=2)
                print(f"Progress: {total} processed")

        else:
            n_not_saved += 1

    print("\nResults:")
    print(f"Total: {total}")
    print(f"correct label matches: {count_label_match} ({100*count_label_match/total:.1f}%)")
    print(f"correct precursor matches: {count_precursor_match} ({100*count_precursor_match/total:.1f}%)")

    print("\nCase Counts:")
    for case_label in CASES:
        count = case_counts.get(case_label, 0)
        print(f"{case_label}: {count} ({100*count/total:.1f}%)")

    if save_dataset:
        print(f"Dataset size: {len(dataset)} ({100*len(dataset)/total:.1f}%)")
        print(f"not saved instances: {n_not_saved} ({100*n_not_saved/total:.1f}%)")

        with open(output_path, 'w') as f:
            json.dump(dataset, f, indent=2)
        print(f"Saved {len(dataset)} training examples to {output_path}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate the reasoning corpus for LLM to perform de novo peptide sequencing")
    parser.add_argument("--input_data_path",
                         default="/disk/raptor-array/mshab008/msbenchmarking/hardness_index_0/filtered_spectra/hardness_index_0.parquet")
    parser.add_argument("--output_path",
                         default="v2/output/sample_llm_reasoning_corpus.json")
    parser.add_argument("--num_examples", type=int, default=10,
                         help="Number of spectra to keep in the output dataset.")
    parser.add_argument("--no_save", action="store_true",
                         help="Process spectra without writing a dataset file (stats-only run).")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    input_data_path = args.input_data_path
    OUTPUT_PATH = args.output_path
    VERBOSE = args.verbose
    NUM_EXAMPLES = args.num_examples if args.num_examples > 0 else None
    SAVE_DATASET = not args.no_save

    process_input_data(input_data_path, OUTPUT_PATH, NUM_EXAMPLES, SAVE_DATASET)