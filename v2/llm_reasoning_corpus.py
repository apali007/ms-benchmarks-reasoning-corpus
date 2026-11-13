import os
import json
import re
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import argparse
from collections import Counter, defaultdict
from datetime import datetime

from constants import (
    AA_MASSES, NEUTRAL_LOSS_MASSES, MASS_PROTON,
    DATASET_METADATA_ROOT, VARIABLE_MODS, SEARCH_OUTPUT_ROOT
)

def tokenize_peptide(modified_peptide):
    # One token per residue, keeping its modification, GM[15.9949]K -> ['G', 'M[15.9949]', 'K']
    return re.findall(r'[A-Z](?:\[[^\]]+\])*', modified_peptide)

def add_modified_residues(aa_masses, tokens):
    #Add a mass entry for each modified residue in this peptide
    masses = dict(aa_masses)
    for token in tokens:
        if '[' in token:
            deltas = re.findall(r'\[([^\]]+)\]', token)
            masses[token] = aa_masses[token[0]] + sum(VARIABLE_MODS[delta][1] for delta in deltas)

    return masses

def residue_parts(residue, fixed_mods):
    # Mass components of a residue: [(label, mass)], the unmodified residue first
    aa = residue[0]
    parts = [(aa, AA_MASSES[aa])]
    if aa in fixed_mods:
        parts.append(("its fixed modification", fixed_mods[aa]))
    for delta in re.findall(r'\[([^\]]+)\]', residue):
        name, mass = VARIABLE_MODS[delta]
        parts.append((name, mass))
    return parts
                            
def residue_text(residue, fixed_mods):
    # Spell out any modifications so the residue mass can be checked against the standard table.
    parts = residue_parts(residue, fixed_mods)
    if len(parts) == 1:
        return residue
    labels = [label for label, _ in parts[1:]]
    mass_text = " + ".join(f"{mass:.5f}" for _, mass in parts)
    return f"{parts[0][0]} with {' and '.join(labels)} ({mass_text} = {sum(m for _, m in parts):.5f} Da)"

def build_ladder(ion_type, matched_ion, label_peptide, precursor_mz, precursor_charge, fixed_mods):
    # Numbers only: one step per observed ion, plus one a precursor-derived terminal residue
    precursor_neutral_mass = precursor_mz * precursor_charge - precursor_charge * MASS_PROTON
    residues = list(reversed(label_peptide)) if ion_type == 'y' else list(label_peptide)
    n = len(residues)
    steps = []
    prev_mass = None
    ptms_explained = set()

    # Position from N-terminus (1-based), so y and b refer to the same residue numbering
    def to_pos(i):
        return n - i + 1 if ion_type == 'y' else i

    for i in range(1, n + 1):
        peak_information = matched_ion.get(f"{ion_type}{i}")

        if peak_information is None:
            # The full length ion is never observed; recover the last residue from the precursor
            if i == n and prev_mass is not None:
                if ion_type == 'y':
                    aa_mass = precursor_neutral_mass - prev_mass + MASS_PROTON
                else:
                    aa_mass = precursor_neutral_mass - NEUTRAL_LOSS_MASSES['H2O'] - prev_mass + MASS_PROTON
                steps.append({
                    "method": "precursor", "ion_type": ion_type, "prev_ion": f"{ion_type}{n-1}",
                    "prev_mass": prev_mass, "precursor_mz": precursor_mz, "precursor_charge": precursor_charge,
                    "neutral_mass": precursor_neutral_mass, "residue_mass": aa_mass,
                    "residue": residues[i-1],
                    "pos": to_pos(i)
                })
    
            prev_mass = None
            continue

        (mz, intensity, charge, loss, tmz) = peak_information
        step = {
            "ion": f"{ion_type}{i}", "mz": float(mz), "theoretical_mz": float(tmz), 
            "intensity": float(intensity), "charge": charge, "loss": loss,
        }
        ion_mass = mz * charge - (charge - 1) * MASS_PROTON
        
        if charge > 1:
            step["singly_charged_mass"] = ion_mass

        if loss is not None:
            ion_mass += NEUTRAL_LOSS_MASSES[loss]
        step['ion_mass'] = ion_mass

        # Name the PTMs on the first observed ion that contains them, as later ions contain them too.
        new_ptms = []
        for j, token in enumerate(residues[:i], start=1):
            if j in ptms_explained:
                continue
            ptms_explained.add(j)
            new_ptms += re.findall(r'\[([^\]]+)\]', token)

        if new_ptms:
            step["ptms"] = [(VARIABLE_MODS[d][0], VARIABLE_MODS[d][1], c) for d, c in Counter(new_ptms).items()]

        if i == 1: 
            if ion_type == 'y':
                aa_mass = ion_mass - NEUTRAL_LOSS_MASSES['H2O'] - MASS_PROTON
            else:
                aa_mass = ion_mass - MASS_PROTON
            step.update(method="first", residue_mass=aa_mass, residue=residues[0], pos=to_pos(1))            
        elif prev_mass is not None:
            step.update(method="diff", prev_ion=f"{ion_type}{i-1}", prev_mass=prev_mass,
            residue_mass=ion_mass - prev_mass, residue=residues[i-1], pos=to_pos(i))
        else:
            #Gap: no previous ion, so no residue named
            step["method"] = None

        steps.append(step)
        prev_mass = ion_mass

    return steps


def render_ladder(steps, explain_precursor, fixed_mods):
    # Text only: every number comes from the steps
    traces = []
    for s in steps:
        if s["method"] == "precursor":
            if explain_precursor:
                precursor_text = (f"The neutral peptide mass is computed from the precursor m/z ({s['precursor_mz']:.5f}) and charge "
                                    f"({s['precursor_charge']}+): {s['precursor_mz']:.5f} x {s['precursor_charge']} - {s['precursor_charge']} x "
                                    f"{MASS_PROTON:.5f} = {s['neutral_mass']:.5f} Da. "
                                    )
            else:
                precursor_text = ""
            if s["ion_type"] == 'y':
                traces.append(
                    f"The remaining residue at the N-terminus comes from the precursor. {precursor_text}"
                    f"The neutral peptide mass ({s['neutral_mass']:.5f} Da) minus {s['prev_ion']} ({s['prev_mass']:.5f} Da) plus one proton "
                    f"({MASS_PROTON:.5f} Da) gives {s['residue_mass']:.5f} Da, which matches the residue mass of {residue_text(s['residue'], fixed_mods)}."
                )
            else:
                traces.append(
                    f"The remaining residue at the C-terminus comes from the precursor. {precursor_text}"
                    f"The neutral peptide mass ({s['neutral_mass']:.5f} Da) minus H2O ({NEUTRAL_LOSS_MASSES['H2O']:.5f} Da), " 
                    f"minus {s['prev_ion']} ({s['prev_mass']:.5f} Da), plus one proton ({MASS_PROTON:.5f} Da) gives {s['residue_mass']:.5f} Da, "
                    f"which matches the residue mass of {residue_text(s['residue'], fixed_mods)}."
                )
                    
            continue

        ion, mz, charge = s["ion"], s["mz"], s["charge"]
        if charge > 1:
            protons_text = f"{MASS_PROTON:.5f}" if charge == 2 else f"{charge-1} x {MASS_PROTON:.5f}"
            line = (f"{ion} is observed at m/z {mz:.5f} as a {charge}+ ion. "
                    f"Converting it to a singly charged ion gives {mz:.5f} x {charge} - {protons_text} = {s['singly_charged_mass']:.5f} Da.")
        else:
            line = f"{ion} is observed at m/z {mz:.5f}."
        
        if s["loss"] is not None:
            line += (f" This peak has lost {s['loss']} ({NEUTRAL_LOSS_MASSES[s['loss']]:.5f} Da), "
                     f"so adding it back gives {s['ion_mass']:.5f} Da.")

        if "ptms" in s:
            parts = []
            for name, mass, count in s["ptms"]:
                if count == 1:
                    parts.append(f"{name} (+{mass:.5f} Da)")
                else:
                    parts.append(f"{count} {name}s ({count} x {mass:.5f} = {count * mass:.5f} Da)")
            line += f" This ion carries {' and '.join(parts)}."
        
        if s["method"] == "first":
            if ion[0] == 'y':
                line += (f" Since y1 contains only the last residue plus H2O and a proton, subtracting H2O "
                f"({NEUTRAL_LOSS_MASSES['H2O']:.5f} Da) and a proton ({MASS_PROTON:.5f} Da) leaves {s['residue_mass']:.5f} Da, "
                f"which matches the residue mass of {residue_text(s['residue'], fixed_mods)}.")
            else:
                line += (f" Since b1 contains only the first residue plus a proton, subtracting the proton "
                f"({MASS_PROTON:.5f} Da) leaves {s['residue_mass']:.5f} Da, which matches the residue mass of {residue_text(s['residue'], fixed_mods)}.")
        elif s["method"] == "diff":
            line += (f" The mass difference between {s['prev_ion']} and {ion} is {s['ion_mass']:.5f} - {s['prev_mass']:.5f} = "
            f"{s['residue_mass']:.5f} Da, which matches the residue mass of {residue_text(s['residue'], fixed_mods)}.")

        traces.append(line)

    return traces

def round_floats(value, digits=5):
    
    # Round every float in the steps for savingl the calculations above keep full precision
    if isinstance(value, float):
        return round(value, digits)
    if isinstance(value, (list, tuple)):
        return [round_floats(v, digits) for v in value]
    if isinstance(value, dict):
        return {k: round_floats(v, digits) for k,v in value.items()}
    
    return value


def ions_and_residues(steps):
    # The record fields the stats already use, derived from the steps
    ions = [{"ion": s["ion"], "mz": round(s["mz"], 5), "theoretical_mz": round(s["theoretical_mz"], 5), 
            "intensity": round(s["intensity"], 5), "charge": s["charge"], "loss": s["loss"]}
            for s in steps if s["method"] != "precursor"]
    residues = [{"pos": s["pos"], "method": s["method"]} for s in steps if s["method"]]

    return ions, residues    

def derive_peptide(matched_y, matched_b, label_peptide, precursor_mz, precursor_charge, fixed_mods):
    y_steps = build_ladder('y', matched_y, label_peptide, precursor_mz, precursor_charge, fixed_mods)
    b_steps = build_ladder('b', matched_b, label_peptide, precursor_mz, precursor_charge, fixed_mods)

    # Explain the precursor conversion in the b ladder only if the y_ladder did not already use it.
    y_used_precursor = any(s['method'] == 'precursor' for s in y_steps)
    y_traces = render_ladder(y_steps, explain_precursor=True, fixed_mods=fixed_mods)
    b_traces = render_ladder(b_steps, explain_precursor=not y_used_precursor, fixed_mods=fixed_mods)
                    

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
        "y_steps": round_floats(y_steps),
        "b_steps": round_floats(b_steps),
    }

    return label_peptide, reasoning_trace, ladder_info
    
def spectrum_to_input_text(mz_list, intensity_list, precursor_mz, precursor_charge, fixed_mods, variable_mods, fragment_tolerance):
    unique_peaks = sorted(set(zip(mz_list, intensity_list)), key=lambda p: p[0])
    # peaks = sorted(zip(mz_list, intensity_list), key=lambda p: p[0])
    mods_text = ", ".join(f"{aa} +{delta:.5f} Da" for aa, delta in fixed_mods.items()) or "none"
    lines = [f"Precursor m/z: {precursor_mz:.5f}, Precursor charge: {precursor_charge}+"]
    lines.append(f"Fixed modifications: {mods_text}")
    lines.append(f"Variable modifications: {', '.join(variable_mods) or 'none'}")
    lines.append(f"Fragment mass tolerance: {fragment_tolerance:g} Da")
    lines += ["Experimental MS/MS spectrum (m/z, intensity):"]
    lines += [f"m/z {mz:.5f}, intensity {intensity:.5f}" for mz, intensity in unique_peaks]
    lines.append("Analyze this spectrum and identify the most likely peptide sequence.")
    return "\n".join(lines)

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

def load_variable_mods(dataset):
        # The metadata JSON has no variable modifications, so read them from the Tide search settings
    params_path = os.path.join(SEARCH_OUTPUT_ROOT, dataset, 'tide-output', 'combined', 'tide-search.params.txt')
    with open(params_path) as f:
        settings = dict(line.strip().split('=', 1) for line in f if '=' in line and not line.startswith('#'))

    sites = {
        'mods-spec': None,
        'nterm-peptide-mods-spec': 'N-terminus',
        'cterm-peptide-mods-spec': 'C-terminus',
        'nterm-protein-mods-spec': 'protein-N-terminus',
        'cterm-protein-mods-spec': 'protein C-terminus',
    }

    variable_mods = []
    for key, site in sites.items():
        for entry in settings.get(key, '').split(','):
            # A leading number marks a variable modification; without one it is static
            m = re.match(r'^\d+([A-Z]+)([+-][\d.]+)$', entry.strip())
            if m is None:
                continue
            residues, delta = m.groups()
            if site is None:
                target = residues
            else:
                target = site if residues == 'X' else f"{residues} at {site}"
            variable_mods.append(f"{target} {float(delta):+.5f} Da")

    return variable_mods

def load_search_params(dataset):

    metadata_path = os.path.join(DATASET_METADATA_ROOT, f"{dataset}.json")
    with open(metadata_path) as f:
        metadata = json.load(f)

    aa_masses = dict(AA_MASSES)
    fixed_mods = {}
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
            fixed_mods[aa] = fixed_mods.get(aa, 0.0) + float(delta)


    return {
        'aa_masses': aa_masses,
        'fragment_tolerance': float(metadata['Fragment_Tolerance']),
        'precursor_window': float(metadata['Precursor_mz_Tolerance']),
        'precursor_window_type': metadata['Precursor_mz_Unit'].strip().lower(),
        'fixed_mods': fixed_mods,
        'variable_mods': load_variable_mods(dataset)
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

            spectra.append({
                "dataset": r.dataset,
                "psmid": r.PSMId,
                "peptide": r.peptide,
                "modified_peptide": r.modified_peptide,
                "precursor_mz": r.precursor_mz,
                "precursor_charge": int(r.charge),
                "mz": r.mz_array.tolist(),
                "intensity": r.intensity_array.tolist(),
            })
            
            n_selected += 1

            if num_examples is not None and n_selected == num_examples:
                print(f"Selected {n_selected}, including {n_modified} modified peptides")
                return spectra

    print(f"Selected {n_selected}, including {n_modified} modified peptides")
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
                candidates.append((loss_rank, charge, error, -intensity, name, loss, peak_idx, tmz))

    # Intact before neutral loss, lower charge first, then closest, then most intense; 
    # each ion and each peak used at most once
    candidates.sort()
    assigned, assigned_peaks = {}, set()
    for loss_rank, charge, error, _, name, loss, peak_idx, tmz in candidates:
        if (name, charge, loss) in assigned or peak_idx in assigned_peaks:
            continue
        assigned[(name, charge, loss)] = (*peaks[peak_idx], tmz)
        assigned_peaks.add(peak_idx)

    # One peak per ladder position: prefer intact over neutral loss, then the lowest charge
    matched_b, matched_y = {}, {}
    for (name, charge, loss), (mz, intensity, tmz) in sorted(assigned.items(), key=lambda kv: (kv[0][2] is not None, kv[0][1])):
        target = matched_b if name.startswith('b') else matched_y
        if name not in target:
            target[name] = (mz, intensity, charge, loss, tmz)

    noise = [peak for idx, peak in enumerate(peaks) if idx not in assigned_peaks]

    return matched_y, matched_b, assigned, noise

def spectrum_stats(modified_peptide, ladder_info):
    tokens = tokenize_peptide(modified_peptide)
    n = len(tokens)
    y_ions, y_residues = ions_and_residues(ladder_info['y_steps'])
    b_ions, b_residues = ions_and_residues(ladder_info['b_steps'])
    n_b, n_y = len(b_ions), len(y_ions)
    supported = {r['pos'] for r in y_residues + b_residues}
    errors = [abs(ion['mz'] - ion['theoretical_mz']) for ion in b_ions + y_ions]

    return {
        'b_coverage': 100 * n_b / (n-1),
        'y_coverage': 100 * n_y / (n-1),
        'ion_coverage': 100 * (n_b + n_y) / (2 * (n-1)),
        'residue_coverage': 100 * len(supported) / n,
        'fully_supported': len(supported) == n,
        'fragment_errors': errors,
        'modified': '[' in modified_peptide,
    }

def print_tables(stats, rows):
    # rows: (label, key in stats) pairs, one table row each
    print("\nIon coverage (observed / expected b and y ions):")
    print(f"{'dataset':10s} {'n':>6s} {'b mean':>8s} {'b median':>9s} {'b complete':>11s} "
            f"{'y mean':>8s} {'y median':>9s} {'y complete':>11s} {'b+y mean':>9s}"
            )
    for label, key in rows:
        s = stats[key]
        b, y = np.array(s['b_coverage']), np.array(s['y_coverage'])
        print(f"{label:10s} {len(b):6d} {b.mean():7.1f}% {np.median(b):8.1f}% {100 * (b == 100).mean():10.1f}% "
                f"{y.mean():7.1f}% {np.median(y):8.1f}% {100 * (y == 100).mean():10.1f}% {np.mean(s['ion_coverage']):8.1f}%"          
                )

    print("\nResidue coverage and fragment mass error:")
    print(f"{'dataset':10s} {'n':>6s} {'residue cov':>12s} {'fully supported':>16s} {'err med':>9s} {'err p95':>9s}")
    for label, key in rows:
        s = stats[key]
        errors = np.array(s['fragment_errors'])
        print(f"{label:10s} {len(s['residue_coverage']):6d} {np.mean(s['residue_coverage']):11.1f}% "
                f"{100 * np.mean(s['fully_supported']):15.1f}% {np.median(errors):9.4f} {np.percentile(errors, 95):9.4f}"
                )


def print_stats(stats):
    datasets = sorted(k for k in stats if isinstance(k, str) and k != 'all')
    print_tables(stats, [(name, name) for name in datasets + ['all']])

    # Compared within a dataset only, because datasets differ in more than their PTMs
    for dataset in datasets:
        if (dataset, 'modified') not in stats:
            continue
        print(f"\n==== {dataset}: unmodified vs modified ====")
        groups = [group for group in ('unmodified', 'modified') if (dataset, group) in stats]
        print_tables(stats, [(group, (dataset, group)) for group in groups])

def add_stats(stats, dataset, spec_stats):
    group = 'modified' if spec_stats['modified'] else 'unmodified'
    for key in (dataset, 'all', (dataset, group)):
        for metric, value in spec_stats.items():
            if metric == 'fragment_errors':
                stats[key][metric].extend(value)
            else:
                stats[key][metric].append(value)


def process_input_data(input_data_path, output_path, num_examples, save_dataset):
    spectra = parse_parquet(input_data_path, num_examples)
    out = open(output_path, 'w') if save_dataset else None
    params_by_dataset = {}
    stats = defaultdict(lambda: defaultdict(list))

    for total, spec in enumerate(spectra, start=1):
        peptide = spec['peptide']
        dataset_name = spec['dataset']
        if dataset_name not in params_by_dataset:
            params_by_dataset[dataset_name] = load_search_params(dataset_name)
        params = params_by_dataset[dataset_name]

        tokens = tokenize_peptide(spec['modified_peptide'])
        aa_masses = add_modified_residues(params['aa_masses'], tokens)

        matched_y, matched_b, annotations, noise = match_ions(
            tokens, spec['mz'], spec['intensity'], 
            aa_masses, params['fragment_tolerance'], spec['precursor_charge'])
        
        final_sequence, reasoning_trace, ladder_info = derive_peptide(
            matched_y,
            matched_b,
            tokens,
            spec['precursor_mz'],
            spec['precursor_charge'],
            params['fixed_mods']
        )

        add_stats(stats, dataset_name, spectrum_stats(spec['modified_peptide'], ladder_info))

        if VERBOSE:
            print(f"\n==== {spec['psmid']} | {spec['modified_peptide']} ====")
            print('\n'.join(reasoning_trace))

        if save_dataset:
            record = {
            'id': spec['psmid'],
            'dataset': dataset_name,
            'spectrum': spectrum_to_input_text(
                spec['mz'], 
                spec['intensity'], 
                spec['precursor_mz'], 
                spec['precursor_charge'],
                params['fixed_mods'],
                params['variable_mods'],
                params['fragment_tolerance']
                ),
            'reasoning_trace': '\n'.join(reasoning_trace),
            'peptide': peptide,
            'modified_peptide': spec['modified_peptide'],
            **ladder_info,
            }
            out.write(json.dumps(record) + '\n')

        if total % 1000 == 0:
            print(f"Progress: {total} processed")

    print(f"\nProcessed {len(spectra)} spectra")
    if save_dataset:
        out.close()
        print(f"Saved {len(spectra)} training examples to {output_path}")

    print_stats(stats)

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
    OUTPUT_PATH = os.path.join(args.output_dir, f"{timestamp}_{n_tag}_llm_reasoning_corpus.jsonl")

    process_input_data(input_data_path, OUTPUT_PATH, NUM_EXAMPLES, SAVE_DATASET)