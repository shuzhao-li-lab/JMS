'''
Based on mass2chem.formula 

1. to find formula for anchor ions
2. to search for any formula (to-do)
'''

from mass2chem.formula import compute_adducts_formulae
from .search import build_centurion_tree, find_best_match_centurion_indexed_list

def get_formula_ions_tree(list_formula_mass, mode='pos'):
    peak_list = []
    for formula, mass in list_formula_mass:
        for ion in compute_adducts_formulae(mass, formula, mode, primary_only=True):
            peak_list.append( {'mz': ion[0], 
                               'neutral_formula': formula, 
                               'ion_relation': ion[1],
                               'neutral_formula_mass': mass} )

    return build_centurion_tree(peak_list)

def search_mz_formula_tree(mz, formula_tree, limit_ppm=5):
    '''
    return the best matched ion (as peak format) in formula_tree.

    formula_tree = get_formula_ions_tree(list_formula_mass, mode='pos')
    '''
    return find_best_match_centurion_indexed_list(mz, formula_tree, limit_ppm)

