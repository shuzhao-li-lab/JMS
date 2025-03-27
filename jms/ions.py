'''
Functions related to isotopes and adductions.
Using mass2chem.

compute_adducts_formulae(304.2038, 'C19H28O3',  mode='pos', primary_only=True)
Out[5]:
[[304.203251, 'M[1+]', 'C19H28O3'],
 [305.21107646677, 'M+H[1+]', 'C19H29O3'],
 [327.19307646677, 'M+Na[1+]', 'C19H28NaO3'],
 [323.22167646677, 'M+H2O+H[1+]', 'C19H31O4']]

'''

from khipu.utils import adduct_search_patterns, adduct_search_patterns_neg, isotope_search_patterns, extended_adducts
from mass2chem.formula import compute_adducts_formulae
from mass2chem.formula import parse_chemformula_dict
import json
import heapq
from scipy.stats import multinomial
import numpy as np
import itertools
from copy import deepcopy

dev_isotopic_patterns = [
    # mass diff, isotopes, (intensity ratio constraint)
    (1.003355, '13C/12C', (0, 0.8)),      # 13C-12C, 12C~99%, 13C ~ 1%
    (0.997035, '15N/14N', (0, 0.2)),     # 15N-14N, 14N ~ 99.64%, 15N ~ 0.36%
    (2.004245, '18O/16O', (0, 0.2)),      # 18O-16O, 16O ~ 99.76, 16O ~ 0.2%
    (1.995796, '34S/32S', (0, 0.4)),      # 32S (95.02%), 33S (0.75%), 34S (4.21%)
    (0.999388, '33S/32S', (0, 0.1)),
    # double isotopes
    (2.00039, 'M(13C),M(15N)', (0, 0.2)),
    (2.999151, 'M(13C),M(34S)', (0, 0.4)),
    # double charged
    (0.5017, '13C/12C, double charged', (0, 0.8)),
    (0.4985, '15N/14N, double charged', (0, 0.2)),
]

#
# -----------------------------------------------------------------------------
#


# todo - we can keep this for now, but this should be reactored 
# alternatively, the Tipu implementation of the isotopologue_generator could be used here.

isotope_information = json.load(open("/Users/mitchjo/Projects/xasari/NIST_isotope_data.json"))

element_to_iso = {}
for iso_dict in isotope_information:
    if iso_dict["Atomic Symbol"] == "D":
        iso_dict["Atomic Symbol"] = "H"
    if iso_dict["Atomic Symbol"] not in element_to_iso:
        element_to_iso[iso_dict["Atomic Symbol"]] = []
    element_to_iso[iso_dict["Atomic Symbol"]].append(iso_dict)

element_to_iso_tuples = {}
for element, isotopes in element_to_iso.items():
    iso_NAPs, iso_masses = [], []
    for iso in isotopes:
        iso_symbol = iso["Atomic Symbol"] + str(iso["Mass Number"])
        mass = iso["Relative Atomic Mass"]
        NAP = iso['Isotopic Composition']
        if NAP and NAP > .005:
            iso_NAPs.append((iso_symbol, NAP, mass))
    iso_NAPs = sorted(iso_NAPs, key=lambda x: -x[1])
    if iso["Atomic Symbol"] == "Sn" or iso["Atomic Symbol"] == "Os":
        iso_NAPs = iso_NAPs[:2]
    if iso_NAPs:
        element_to_iso_tuples[element] = iso_NAPs

multi_cache = {}
element_prob_vectors = {}
def calc_prob(iso_coord):
    sub_coord_probs = []
    for element, sub_coord in iso_coord.items():
        key_tuple = (element, tuple(sub_coord))
        if key_tuple in multi_cache:
            sub_coord_probs.append(multi_cache[key_tuple])
        else:
            key_sum = (element, sum(sub_coord))
            if element not in element_prob_vectors:
                element_prob_vectors[element] = np.array([x[1] for x in element_to_iso_tuples[element]])
            if key_sum not in multi_cache:
                multi_cache[key_sum] = multinomial(key_sum[1], element_prob_vectors[element])
            multi_cache[key_tuple] = multi_cache[key_sum].pmf(sub_coord)
            sub_coord_probs.append(multi_cache[key_tuple])
    return np.product(sub_coord_probs)

mass_cache = {}
element_mass_vectors = {}
def calc_mass(iso_coord):
    sub_coord_masses = []
    for element, sub_coord in iso_coord.items():
        key_tuple = (element, tuple(sub_coord))
        if key_tuple in mass_cache:
            sub_coord_masses.append(mass_cache[key_tuple])
        else:
            if element not in element_mass_vectors:
                element_mass_vectors[element] = np.array([x[2] for x in element_to_iso_tuples[element]])
            mass_cache[key_tuple] = np.dot(sub_coord, element_mass_vectors[element])
        sub_coord_masses.append(mass_cache[key_tuple])
    return np.sum(sub_coord_masses)

neighbor_cache = {}
delta_cache = {}
def gen_neighbors(iso_coord):
    new_coords = []
    for _, sub_coord in iso_coord.items():
        key = tuple(sub_coord)
        if key not in neighbor_cache: 
            neighbor_cache[key] = [sub_coord]
            if len(sub_coord) > 1:
                if len(sub_coord) not in delta_cache:
                    delta_cache[len(sub_coord)] = np.array([x for x in itertools.permutations([-1, 1] + [0 for _ in range(len(sub_coord)-2)])], dtype=np.int16)           
                for delta in delta_cache[len(sub_coord)]:
                    new_coord = delta + sub_coord
                    if not np.any(new_coord < 0):
                        neighbor_cache[key].append(new_coord)
        new_coords.append(neighbor_cache[key])
    new_iso_coords = []
    for coord_mask in np.ndindex(tuple([len(x) for x in new_coords])):
        if not np.all(coord_mask == 0):
            new_iso_coord = {}
            for i, (index, k) in enumerate(zip(coord_mask, iso_coord.keys())):
                new_iso_coord[k] = new_coords[i][index]
            new_iso_coords.append(new_iso_coord)
    return new_iso_coords

def hash_isotopologue(isotopologue):
    return "|".join([k + ','.join([str(x) for x in v]) for k,v in isotopologue.items()])

delta_substrings = {}
isotopes_vectors = {}
def delta_info(isotopologue, reference):
    delta_string = []
    for element in isotopologue.keys():
        iso_sub_coord, ref_sub_coord = isotopologue[element], reference[element]
        key = (element, tuple(iso_sub_coord), tuple(ref_sub_coord))
        if key not in delta_substrings:
            substring = []
            delta_values = iso_sub_coord - ref_sub_coord
            if element not in isotopes_vectors:
                isotopes_vectors[element] = [x[0] for x in element_to_iso_tuples[element]]
            for delta, iso_name in zip(delta_values, isotopes_vectors[element]):
                if delta > 1:
                    substring.append(str(delta) + iso_name)
                elif delta == 1:
                    substring.append(iso_name)
                else:
                    pass
                delta_substrings[key] = ",".join(substring)
        delta_string.append(delta_substrings[key])
    delta_string = '(' + ",".join([x for x in delta_string if x]) + ')'
    if delta_string != '()':
        return delta_string, calc_mass(isotopologue) - calc_mass(reference)
    else:
        return '', calc_mass(isotopologue) - calc_mass(reference)


def gen_isotopologues(formula, NAP_cutoff=0):
    iso_coord = {}
    for k, v in parse_chemformula_dict(formula).items():
        iso_coord[k] = [0 for _ in element_to_iso_tuples[k]]
        iso_coord[k][0] = v
        iso_coord[k] = np.array(iso_coord[k], dtype=np.int16)
    iso_heap = []
    isos = [(-1 * calc_prob(x), x) for x in [iso_coord] + gen_neighbors(iso_coord)]
    isos = sorted(isos, key=lambda x: x[0])
    reference = isos[0][1]
    used = set()
    for NAP_prob, iso_coord in isos:
        if -1 * NAP_prob > NAP_cutoff:
            iso_hash = hash_isotopologue(iso_coord)
            if iso_hash not in used:
                heapq.heappush(iso_heap, (NAP_prob, iso_coord))
                used.add(hash_isotopologue(iso_coord))
    while iso_heap:
        NAP, most_abundant = heapq.heappop(iso_heap)
        yield -1 * NAP, most_abundant, calc_mass(most_abundant), delta_info(most_abundant, reference)
        for x in gen_neighbors(most_abundant):
            hash = hash_isotopologue(x)
            if hash not in used:
                NAP_prob = calc_prob(x)
                if NAP_prob > NAP_cutoff:
                    heapq.heappush(iso_heap, (-1 * NAP_prob, x))
                    used.add(hash)

def generate_ion_signature(mw, neutral_formula,  mode='pos', primary_only=True):
    '''
    Extend mass2chem.formula.compute_adducts_formulae by C13.
    Note - Resulting chemical formula is not computable.
    '''
    adducts = compute_adducts_formulae(mw, neutral_formula,  mode, primary_only)
    C13 = []
    for A in adducts:
        C13.append([
            A[0]+1.003355, A[1]+',C13', A[2]+',(C13)'
        ])
    return adducts + C13

def generate_ion_signature_comprehensive(mw, neutral_formula, mode='pos', primary_only=True, C13_only=False, NAP_cutoff=.01):
    '''
    todo - write docstring
    '''
    adducts = compute_adducts_formulae(mw, neutral_formula, mode, primary_only)
    if C13_only:
        C13 = []
        for A in adducts:
            A.append(0)
            C13.append([
                A[0]+1.003355, A[1]+',C13', A[2]+',(C13)', 1
            ])
        return adducts + C13
    else:
        isotopologues = []
        for i, (NAP, _, _, (delta_string, delta_mass)) in enumerate(gen_isotopologues(neutral_formula)):
            if NAP > NAP_cutoff:
                for A in adducts:
                    if delta_string:
                        isotopologues.append([
                            A[0] + delta_mass,
                            A[1] + "," + delta_string[1:-1] + ";" + str(i),
                            A[2] + "," + delta_string, 
                            i,
                        ])
                    else:
                        isotopologues.append([
                            A[0] + delta_mass,
                            A[1] + ";" + str(i),
                            A[2],
                            i,
                        ])
            else:
                return isotopologues
        isotopologues = []
        for A in adducts:
            isotopologues.append([
                A[0],
                A[1] + ";0",
                A[2],
                0,
            ])
        return isotopologues
