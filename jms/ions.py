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

