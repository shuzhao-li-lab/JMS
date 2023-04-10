'''
Based on mass2chem.formula 

1. to find formula for anchor ions
2. to search for any formula (to-do)
'''

from mass2chem.formula import compute_adducts_formulae, \
                                    parse_chemformula_dict, \
                                        add_formula_dict, \
                                            dict_to_hill_formula

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

def adjust_charge_in_formula(charged_formula:str,charge:int)->str: 
    """adjust charged formula to neutral formula

    Parameters
    ----------
    charged_formula : str
        charged formula
    charge : int
        charge

    Returns
    -------
    str
        neutral formula

    Examples
    --------
    >>> adjust_charge_in_formula('C17H31O2',-1)
    'C17H32O2'
    
    """ 
    if charge == 0:
        result = charged_formula
    else:
        formula_dict = parse_chemformula_dict(charged_formula)
        hydrogen_dict = {'H':-charge}
        hypothetical_neutral_formula_dict = add_formula_dict(formula_dict, hydrogen_dict)
        if hypothetical_neutral_formula_dict:
            result = dict_to_hill_formula(hypothetical_neutral_formula_dict)
        else: 
            if charge < 0:
                result = None
            else:  # if charge >0, e.g., Fe, charge is 2+; it should still be Fe as formula
                result = charged_formula # here dealing with Metal etc.
    return result
