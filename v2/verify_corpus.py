import argparse
import json
import re
from collections import Counter, defaultdict

import numpy as np
import spectrum_utils.spectrum as sus
from pyteomics import mass

# Neutral loss from pyteomics, independent of constants.py; spectrum_utils wants them as negative masses
NEUTRAL_LOSSES = {'H2O': -mass.calculate_mass(formula='H2O'), 'NH3': -mass.calculate_mass(formula='NH3')}

def parse_spectrum(text):
    # Read precursor, tolerance and peaks back from the spectrum text the model sees
    precursor = re.search(r'Precursor m/z: ([\d.]+), Precursor charge: (\d+)\+', text)
    tolerance = re.search(r'Fragment mass tolerance: ([\d.]+) Da', text)
    fixed = re.search(r'Fixed modifications: (.*)', text).group(1)
    peaks = [(float(mz), float(intensity))
             for mz, intensity in re.findall(r'm/z ([\d.]+), intensity ([\d.]+)', text)]
    return {
        "precursor_mz": float(precursor.group(1)),
        "precursor_charge": int(precursor.group(2)),
        "tolerance": float(tolerance.group(1)),
        "fixed_mods": re.findall(r'([A-Z]) \+([\d.]+) Da', fixed),
        "peaks": peaks,
    }

def to_proforma(modified_peptide, fixed_mods):
    # KEM[15.9949]QPTHPIR with C +57.02146 -> <[+57.02146]@C>KEM[+15.9949]QPTHPIR
    fixed = ''.join(f"<[+{delta}]@{aa}>" for aa, delta in fixed_mods)
    return fixed + re.sub(r'\[([\d.]+)\]', r'[+\1]', modified_peptide)

def annotate(record, spectrum):
    # spectrum_utils labels every peak it can explain: {(m/z, intensity): {(ion, charge, loss): theoretical ions m/z}}
    # spectrum_utils sorts peaks by m/z, so sort them the same way to keep labels aligned with peaks
    peaks = sorted(spectrum['peaks'])
    mz = np.array([p[0] for p in peaks])
    intensity = np.array([p[1] for p in peaks], dtype=np.float32)
    spec = sus.MsmsSpectrum(record['id'], spectrum['precursor_mz'], spectrum['precursor_charge'], mz, intensity)
    spec.annotate_proforma(to_proforma(record['modified_peptide'], spectrum['fixed_mods']),
                           spectrum['tolerance'], 'Da', ion_types='by', 
                           max_ion_charge=spectrum['precursor_charge'], neutral_losses=NEUTRAL_LOSSES)
    labels = defaultdict(dict)
    for peak, interpretation in zip(peaks, spec.annotation):
        for f in interpretation.fragment_annotations:
            loss = f.neutral_loss[1:] if f.neutral_loss else None
            labels[peak][(f.ion_type, f.charge, loss)] = peak[0] - f.mz_delta[0]
    return labels

def check_ions(record, spectrum, fail):
    labels = annotate(record, spectrum)
    peaks = set(spectrum['peaks'])
    for s in record['y_steps'] + record['b_steps']:
        if s['method'] == 'precursor':        
            continue
        key = (s['ion'], s['charge'], s['loss'])
        # 1. the peak is really in the spectrum the model sees
        fail('peak_in_spectrum', (s['mz'], s['intensity']) not in peaks, s)
        # 2. within the fragment tolerance
        fail('within_tolerance', abs(s['mz'] - s['theoretical_mz']) > spectrum['tolerance'] + 1e-5, s)
        # 3. spectrum_utils labels this peak as the same ion, charge and loss
        found = labels.get((s['mz'], s['intensity']), {}).get(key)
        fail('spectrum_utils_same_ion', found is None, s)
        # 4. and computes the same theoretical m/z
        if found is not None:
            fail('theoretical_mz', abs(found - s['theoretical_mz']) > 1e-4, s)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verify the ions and calculations in a reasoning corpus")
    parser.add_argument("corpus_path")
    args = parser.parse_args()

    checked = defaultdict(Counter) # dataset -> check -> number checked
    failed = defaultdict(Counter) # dataset -> check -> number failed

    n = 0
    with open(args.corpus_path) as f:
        for line in f:
            record = json.loads(line)
            spectrum = parse_spectrum(record['spectrum'])
            dataset = record['dataset']

            def fail(check, is_failure, step):
                checked[dataset][check] += 1
                if is_failure:
                    failed[dataset][check] += 1
                    if sum(failed[dataset].values()) <= 5:
                        print(f"FAIL {check} {record['id']} {record['modified_peptide']} {step}")

            check_ions(record, spectrum, fail)
            n += 1

    print(f"Read {n} records from {args.corpus_path}")
    for dataset in sorted(checked):
        print(dataset)
        for check, total in checked[dataset].items():
            print(f" {check:24s} {total - failed[dataset][check]:>8} / {total:<8} passed")
