'''
Convenience functions to construct empirical compounds via khipu;
and retrieve, filter and search empCpds.
'''

import json
import numpy as np
from .search import find_all_matches_centurion_indexed_list
from khipu.epdsConstructor import epdsConstructor
from khipu.utils import adduct_search_patterns, \
                            adduct_search_patterns_neg, \
                                isotope_search_patterns, \
                                    extended_adducts


def get_khipu_epds_from_list_peaks(
                    list_peaks, 
                    mode, 
                    isotope_search_patterns,
                    adduct_patterns,
                    extended_adducts,
                    mz_tolerance_ppm,
                    rt_tolerance):
    '''
    Wrapper function to build empirical compounds from a list of JSON features/peaks.

    Parameters
    ----------
    list_peaks : list of features, [{'mz': 133.097023, 'rtime': 654, 'height': 14388.0, 'id': 555}, ...]
    mode : ionization mode, 'pos' or 'neg'.
    isotope_search_patterns : E.g. [ (1.003355, '13C/12C', (0, 0.8)), (2.00671, '13C/12C*2', (0, 0.8)),]
    adduct_patterns : e.g. [ (21.9820, 'Na/H'), (41.026549, 'Acetonitrile')]
    extended_adducts : adducts used for 2nd round search, e.g. [(117.02655, '-NH3'),
                            (17.02655, 'NH3'), (-18.0106, '-H2O'), ...]
    mz_tolerance_ppm : ppm tolerance in examining m/z patterns.
    rt_tolerance : tolerance of retention time, arbitrary unit but consistent with list_peaks.

    Returns
    -------
    dict_empCpds : {interim_id: empCpd, ...}
    '''
    ECCON = epdsConstructor(list_peaks, mode=mode)
    return ECCON.peaks_to_epdDict(
                    isotope_search_patterns,
                    adduct_patterns,
                    extended_adducts,
                    mz_tolerance_ppm,
                    rt_tolerance,
    ) 

def load_epds_from_json(file):
    '''Read JSON annotation dictionary from asari/khipu.
    Returns list of empCpds.
    '''
    epds = json.load(open(file))
    return list(epds.values())


def filter_epds(list_epds, neutral_formula_mass=True, multiple_ions=True, C13=True):
    '''Filter list of empCpds by neutral_formula_mass or by multiple_ions.
    list_epds : [{'interim_id': 'kp100_128.0951', 'neutral_formula_mass': 128.09508427175538, 'neutral_formula': None, 
        'Database_referred': [], 'identity': [], 'MS1_pseudo_Spectra': [
        {'apex': 653, 'peak_area': 8047955, 'height': 988673, 'left_base': 642, 'right_base': 659, 
        'goodness_fitting': 0.9397139863801345, 'cSelectivity': 0.3497536945812808, 
        'parent_masstrack_id': 880, 'mz': 130.10569763183594, 'snr': 11, 'id_number': 'F1444', 
        'rtime': 125.24219400000001, 'rtime_left_base': 123.24607800000001, 'rtime_right_base': 126.27832200000002, 
        'representative_intensity': 8047955, 'id': 'F1444', 'isotope': '13C/12C', 
        'modification': 'M+H+', 'ion_relation': '13C/12C,M+H+', 'parent_epd_id': 'kp100_128.0951'}, 
        {'apex': 654, 'peak_area': 84934827, 'height': 11980759, 'left_base': 648, 'right_base': 659, 
        'goodness_fitting': 0.993842814809182, 'cSelectivity': 0.43617021276595747, 
        'parent_masstrack_id': 858, 'mz': 129.10237884521484, 'snr': 3, 'id_number': 'F1261', 
        'rtime': 125.41953000000001, 'rtime_left_base': 124.33713, 'rtime_right_base': 126.27832200000002, 
        'representative_intensity': 84934827, 'id': 'F1261', 'isotope': 'M0', 
        'modification': 'M+H+', 'ion_relation': 'M0,M+H+', 'parent_epd_id': 'kp100_128.0951'}], 
        'MS2_Spectra': []}, 
        ...]
    '''
    result = list_epds
    if neutral_formula_mass:
        result = [x for x in result if x['neutral_formula_mass']]
    if multiple_ions:
        result = [x for x in result if len(x['MS1_pseudo_Spectra'])>1]
    if C13:
        result = [x for x in result if check_13C_M1(x)]
    return result

def check_13C_M1(empCpd):
    '''Check presence of pair of M1 13C ion and the M0 12C ion.
    '''
    ions = [x['isotope'] for x in empCpd['MS1_pseudo_Spectra'] if 'isotope' in x]
    if '13C/12C' in ions and 'M0' in ions:
        return True
    else:
        return False

def get_neutrals(list_epds):
    neutrals = []
    for v in list_epds:
        p = {}
        p['id'] = v['interim_id']
        p['mz'] = v['neutral_formula_mass']
        p['rtime'] = np.mean([x['rtime'] for x in v['MS1_pseudo_Spectra']])
        neutrals.append(p)
    return neutrals

def get_match(cpds, mztree, ppm=5):
    '''Find matches of a list of cpds in mztree.
    cpds : [{'id': 'C00025', 'mw': 147.0532, 'name': 'L-Glutamate'}, ...]
    '''
    match = []
    for x in cpds:
        if 'mw' in x and x['mw']:
            _m = find_all_matches_centurion_indexed_list(x['mw'], mztree, ppm)
            if _m:
                match.append( (x, _m) )
    return match
