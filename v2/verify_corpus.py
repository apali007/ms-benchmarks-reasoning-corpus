import argparse
import json
import re
from collections import Counter, defaultdict
import pyopenms as oms
from pyopenms.Constants import PROTON_MASS_U

# masses from OpenMS, independent of constants.py
LOSSES = {'H2O': oms.EmpiricalFormula('H2O').getMonoWeight(), 'NH3': oms.EmpiricalFormula('NH3').getMonoWeight()}

PROTON = PROTON_MASS_U
# steps are saved to 5 decimal places (each off up to 0.5e-5); a few combined, times the charge, stay under 3e-5
ROUNDING = 3e-5

# Amino acid masses; Internal = the residues inside a chain (no water)
RESIDUES = oms.ResidueDB()

# PTM names used in the corpus -> (Unimod name, residue, where); OpenMS looks up their exact masses in Unimod
SITE = oms.ResidueModification.TermSpecificity
UNIMOD_NAMES = {'oxidation': ('Oxidation', 'M', SITE.ANYWHERE), 'acetylation': ('Acetyl', '', SITE.N_TERM)}
UNIMOD_MASSES = {name: oms.ModificationsDB().getModification(*unimod).getDiffMonoMass()
                 for name, unimod in UNIMOD_NAMES.items()}

# OpenMS theoretical b/y ions, set up once; it adds H2O/NH3 losses only for fragments with a residue that can lose them
GENERATOR = oms.TheoreticalSpectrumGenerator()
params = GENERATOR.getParameters()
for option in ('add_b_ions', 'add_y_ions', 'add_first_prefix_ion', 'add_losses', 'add_metainfo'):
    params.setValue(option, 'true')
GENERATOR.setParameters(params)
# corpus loss -> OpenMS ion name part, e.g. y4 with H2O loss at 2+ is 'y4-H2O1++'
OPENMS_LOSSES = {None: '', 'H2O': '-H2O1', 'NH3': '-H3N1'}

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
    variable = re.search(r'Variable modifications: (.*)', text).group(1)
    peaks = [(float(mz), float(intensity))
             for mz, intensity in re.findall(r'm/z ([\d.]+), intensity ([\d.]+)', text)]
    return {
        "precursor_mz": float(precursor.group(1)),
        "precursor_charge": int(precursor.group(2)),
        "tolerance": float(tolerance.group(1)),
        "precursor_tolerance": (float(precursor_tol.group(1)), precursor_tol.group(2)),
        "fixed_mods": re.findall(r'([A-Z]) \+([\d.]+) Da', fixed),
        "peaks": peaks,
        "nterm_mods": [float(d) for d in re.findall(r'N-terminus \+([\d.]+) Da', variable)],
    }

def to_openms(tokens, spectrum):
    # ['L[42.0106]', 'P', 'M[15.9949]', 'C'] -> .[+42.0106]LPM[+15.9949]C[+57.02146]
    # the corpus writes an N-terminal mod on the first residue; OpenMS wants it in front of the sequence
    nterm, residues = '', []
    for pos, token in enumerate(tokens):
        aa, deltas = token[0], re.findall(r'\[([^\]]+)\]', token)
        if pos == 0:
            is_nterm = [any(abs(float(d) - m) < 1e-3 for m in spectrum['nterm_mods']) for d in deltas]
            nterm = ''.join(f"[+{d}]" for d, n in zip(deltas, is_nterm) if n)
            deltas = [d for d, n in zip(deltas, is_nterm) if not n]
        deltas += [delta for fixed_aa, delta in spectrum['fixed_mods'] if fixed_aa == aa]
        residues.append(aa + ''.join(f"[+{d}]" for d in deltas))
    return ('.' + nterm if nterm else '') + ''.join(residues)

def theoretical_ions(tokens, spectrum):
    # OpenMS theoretical ions of the label peptide, charges 1..precursor charge: {'y4++': m/z, 'b2-H2O1+': m/z, ...}
    ions = oms.MSSpectrum()
    GENERATOR.getSpectrum(ions, oms.AASequence.fromString(to_openms(tokens, spectrum)), 1, spectrum['precursor_charge'])
    names = ions.getStringDataArrays()[0]
    return {name.decode(): peak.getMZ() for peak, name in zip(ions, names)}

def check_ions(record, spectrum, tokens, fail):
    theoretical = theoretical_ions(tokens, spectrum)
    peaks = set(spectrum['peaks'])
    tol = spectrum['tolerance']
    for s in record['y_steps'] + record['b_steps']:
        if s['method'] == 'precursor':
            continue
        # 1. the peak is really in the spectrum the model sees
        fail('peak_in_spectrum', (s['mz'], s['intensity']) not in peaks, s)
        # 2. within the fragment tolerance
        fail('within_tolerance', abs(s['mz'] - s['theoretical_mz']) > tol + 1e-5, s)
        # 3. OpenMS has the same ion (type, number, loss, charge) and the peak is within tolerance of it
        found = theoretical.get(f"{s['ion']}{OPENMS_LOSSES[s['loss']]}{'+' * s['charge']}")
        fail('openms_same_ion', found is None or abs(s['mz'] - found) > tol + 1e-5, s)
        # 4. and computes the same theoretical m/z
        if found is not None:
            fail('theoretical_mz', abs(found - s['theoretical_mz']) > 1e-4, s)


def residue_true_mass(token, fixed_mods):
    # Residue mass from OpenMS plus its fixed and variable modifications
    aa = token[0]
    deltas = [float(d) for d in re.findall(r'\[([^\]]+)]', token)]
    return RESIDUES.getResidue(aa).getMonoWeight(oms.Residue.ResidueType.Internal) + float(dict(fixed_mods).get(aa, 0.0)) + sum(deltas)

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
                fail('singly_charged_mass', 'singly_charged_mass' not in s or abs(singly - s['singly_charged_mass']) > ROUNDING, s)
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
            # The residue named is the label residue at pos, and its mass agrees with OpenMs
            pos = n - i + 1 if ion_type == 'y' else i
            fail('residue_label', s['pos'] != pos or s['residue'] != tokens[pos - 1], s)
            true_mass = residue_true_mass(s['residue'], spectrum['fixed_mods'])
            fail('residue_vs_openms', abs(s['residue_mass'] - true_mass) > allowed + ROUNDING, s)


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
        # an ion nmber longer than the peptide is already a failure; only look at real positions
        new &= set(range(1, n+1))
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

            # tokens are split into residues, keeping each mod with its residue
            tokens = re.findall(r'[A-Z](?:\[[^\]]+\])*', record['modified_peptide'])
            check_ions(record, spectrum, tokens, fail)
            for ion_type in ('y', 'b'):
                check_steps(record[f'{ion_type}_steps'], ion_type, tokens, spectrum, fail)
                check_ptms(record[f'{ion_type}_steps'], ion_type, tokens, fail)

            n += 1

    print(f"Read {n} records from {args.corpus_path}")
    for dataset in sorted(checked):
        print(dataset)
        for check, total in checked[dataset].items():
            print(f" {check:24s} {total - failed[dataset][check]:>8} / {total:<8} passed")
