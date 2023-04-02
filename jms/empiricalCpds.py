'''
Construction of empirical compounds via khipu 


Will add mummichog functions on userData and empCpds here -



'''

import json
import numpy as np
from .search import find_all_matches_centurion_indexed_list












def load_epds_from_json(file):
    '''Read JSON annotation dictionary from asari/khipu.
    Returns list of empCpds.
    '''
    epds = json.load(open(file))
    return list(epds.values())

def filter_epds(list_epds, neutral_formula_mass=True, multiple_ions=True):
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
    return result


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



