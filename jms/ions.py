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

from khipu.utils import adduct_search_patterns, \
                            adduct_search_patterns_neg, \
                                isotope_search_patterns, \
                                    extended_adducts


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

from mass2chem.formula import compute_adducts_formulae

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

