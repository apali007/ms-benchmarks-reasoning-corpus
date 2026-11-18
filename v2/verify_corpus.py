import argparse
import json
import re
from collections import Counter, defaultdict
import numpy as np
import spectrum_utils.proforma as proforma
import spectrum_utils.spectrum as sus
from pyteomics import mass

# Neutral loss from pyteomics, independent of constants.py; spectrum_utils wants them as negative masses
NEUTRAL_LOSSES = {'H2O': -mass.calculate_mass(formula='H2O'), 'NH3': -mass.calculate_mass(formula='NH3')}

PROTON = mass.nist_mass['H+'][0][0]
LOSSES = {name: -m for name, m in NEUTRAL_LOSSES.items()}
# steps are saved to 5 decimal places (each off up to 0.5e-5); a few combined, times the charge, stay under 3e-5
ROUNDING = 3e-5

# PTM names used in the corpus -> Unimod names; spectrum_utils looks up their exact masses in Unimod
UNIMOD_NAMES = {'oxidation': 'Oxidation', 'acetylation': 'Acetyl'}
UNIMOD_MASSES = {name: proforma.parse(f"A[{unimod}]")[0].modifications[0].mass
                 for name, unimod in UNIMOD_NAMES.items()}

def precursor_error(neutral_mass, tolerance):
    # Largest error in the neutral mass: ppm scales with the mass, Da is fixed
    value, unit = tolerance
    return neutral_mass * value * 1e-6 if unit == 'ppm' else value

def parse_spectrum(text):
    # Read precursor, tolerance and peaks back from the spectrum text the model sees
    precursor = re.search(r'Precursor m/z: ([\d.]+), Precursor charge: (\d+)\+', text)
    tolerance = re.search(r'Fragment mass tolerance: ([\d.]+) Da', text)
    precursor_tol = re.search(r'Precursor mass tolerance: ([\d.]+) (ppm|Da)', text)
    fixed = re.search(r'Fixed modifications: (.*)', text).group(1)
    peaks = [(float(mz), float(intensity))
             for mz, intensity in re.findall(r'm/z ([\d.]+), intensity ([\d.]+)', text)]
    return {
        "precursor_mz": float(precursor.group(1)),
        "precursor_charge": int(precursor.group(2)),
        "tolerance": float(tolerance.group(1)),
        "precursor_tolerance": (float(precursor_tol.group(1)), precursor_tol.group(2)),
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


def residue_true_mass(token, fixed_mods):
    # Residue mass from pyteomics plus its fixed and variable modifications
    aa = token[0]
    deltas = [float(d) for d in re.findall(r'\[([^\]]+)]', token)]
    return mass.std_aa_mass[aa] +float(dict(fixed_mods).get(aa, 0.0)) + sum(deltas)

def check_steps(steps, ion_type, tokens, spectrum, fail):
    # Recompute every number in the steps from the step itself
    n = len(tokens)
    tol = spectrum['tolerance']
    prev = None # previous ion step

    for s in steps:
        if s['method'] == 'precursor':
            z = spectrum['precursor_charge']
            neutral = spectrum['precursor_mz'] * z - z * PROTON
            fail('neutral_mass', abs(neutral - s['neutral_mass']) > ROUNDING, s)
            if ion_type == 'y':
                residue_mass = s['neutral_mass'] - s['prev_mass'] + PROTON
            else:
                residue_mass = s['neutral_mass'] - LOSSES['H2O'] - s['prev_mass'] + PROTON
            fail('residue_mass', abs(residue_mass - s['residue_mass']) > ROUNDING, s)
            fail('prev_ion', prev is None or s['prev_ion'] != prev['ion'] or s['prev_mass'] != prev['ion_mass']
                 or s['prev_ion'] != f"{ion_type}{n-1}", s)
            i = n
            # previous ion's m/z error times its charge, plus the precursor mass error (ppm of the neutral mass or Da)
            allowed = tol * prev['charge'] + precursor_error(s['neutral_mass'], spectrum['precursor_tolerance']) if prev else 0
        else:
            i = int(s['ion'][1:])
            z = s['charge']
            singly = s['mz'] * z - (z-1) * PROTON
            if z > 1:
                fail('singly_charged_mass', abs(singly - s['singly_charged_mass']) > ROUNDING, s)
            ion_mass = singly + (LOSSES[s['loss']] if s['loss'] else 0.0)
            fail('ion_mass', abs(ion_mass - s['ion_mass']) > ROUNDING, s)

            if s['method'] == 'first':
                residue_mass = s['ion_mass'] - (LOSSES['H2O'] + PROTON if ion_type == 'y' else PROTON)
                fail('residue_mass', abs(residue_mass - s['residue_mass']) > ROUNDING, s)
                # one ion: its m/z error (tol) becomes tol * z in mass
                allowed = tol * z
            elif s['method'] == 'diff':
                fail('prev_ion', prev is None or s['prev_ion'] != prev['ion'] or s['prev_mass'] != prev['ion_mass']
                     or s['prev_ion'] != f"{ion_type}{i-1}", s)
                fail('residue_mass', abs(s['ion_mass'] - s['prev_mass'] - s['residue_mass']) >ROUNDING, s)
                # two ions: both m/z errors, each times its charge, add up
                allowed = tol * (z + prev['charge']) if prev else 0
            prev = s

        if s['method']:
            # The residue named is the label residue at pos, and its mass agrees with pyteomics
            pos = n - i + 1 if ion_type == 'y' else i
            fail('residue_label', s['pos'] != pos or s['residue'] != tokens[pos - 1], s)
            true_mass = residue_true_mass(s['residue'], spectrum['fixed_mods'])
            fail('residue_vs_pyteomics', abs(s['residue_mass'] - true_mass) > allowed + ROUNDING, s)


def check_ptms(steps, ion_type, tokens, fail):
    # Each PTM is named, with its Unimod mass, on the first observed ion that covers its residue
    n = len(tokens)
    covered = set() # positions covered by earlier ions this ladder
    for s in steps:
        if s['method'] == 'precursor':
            continue
        i = int(s['ion'][1:])
        # b_i covers positions 1..i, y_i covers the last i positions
        positions = set(range(1, i+1)) if ion_type == 'b' else set(range(n - i + 1, n + 1))
        # residues this ion is first to cover
        new = positions - covered
        covered |= positions
        # mods expected on this ion: the [delta]s on the newly covered residues, e.g. {'15.9949': 1}
        expected = Counter(d for p in new for d in re.findall(r'\[([^\]]+)\]', tokens[p-1]))

        found = Counter()
        for name, ptm_mass, count in s.get('ptms', []):
            unimod = UNIMOD_MASSES.get(name)
            fail('ptm_mass_unimod', unimod is None or abs(ptm_mass - unimod) > ROUNDING, s)
            # which [delta] in the peptide this named PTM is
            delta = next((d for d in expected if unimod is not None and abs(float(d) - unimod) < 1e-4), None)
            found[delta] += count
        if expected or found:
            fail('ptm_placement', found != expected, s)

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
            # tokens are split into residues, keeping each mod with its residue
            tokens = re.findall(r'[A-Z](?:\[[^\]]+\])*', record['modified_peptide'])
            for ion_type in ('y', 'b'):
                check_steps(record[f'{ion_type}_steps'], ion_type, tokens, spectrum, fail)
                check_ptms(record[f'{ion_type}_steps'], ion_type, tokens, fail)

            n += 1

    print(f"Read {n} records from {args.corpus_path}")
    for dataset in sorted(checked):
        print(dataset)
        for check, total in checked[dataset].items():
            print(f" {check:24s} {total - failed[dataset][check]:>8} / {total:<8} passed")
