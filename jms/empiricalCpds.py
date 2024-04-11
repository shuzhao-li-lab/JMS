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



#
# ------------------------------------------------------------------------------------------
# functions to summarize khipus
#


def get_feature_of_max_intensity(featureList):
    '''
    To get feature of max intensity, and avoid errors by sorting.
    e.g. sorted(M0, reverse=True)[0][1] leas to `TypeError: '<' not supported between instances of 'dict' and 'dict'`
    Use np.argmax here, which is okay with ties.
    '''
    ints = [f['representative_intensity'] for f in featureList]
    idx = np.argmax(ints)
    return featureList[idx]

def get_M0(MS1_pseudo_Spectra):
    '''returns M0 feature with highest representative_intensity.
    Without verifying which ion form.'''
    M0 = [f for f in MS1_pseudo_Spectra if f['isotope']=='M0']
    if M0:
        return get_feature_of_max_intensity(M0)
    else:
        return []
    
def get_M1(MS1_pseudo_Spectra):
    '''returns M+1 feature with highest representative_intensity.
    Without verifying which ion form.'''
    M = [f for f in 
          MS1_pseudo_Spectra if f['isotope']=='13C/12C']
    if M:
        return get_feature_of_max_intensity(M)
    else:
        return []
     
def get_highest_13C(MS1_pseudo_Spectra):
    '''returns 13C labeled feature with highest representative_intensity.
    Without verifying which ion form. Because the label goes with sepecific atoms depending on pathway.
    '''
    M = [f for f in 
          MS1_pseudo_Spectra if '13C/12C' in f['isotope']]
    if M:
        return get_feature_of_max_intensity(M)
    else:
        return []

    
def filter_khipus(list_empCpds, natural_ratio_limit=0.5):
    '''
    returns 
    isopair_empCpds = with good natural 13C ratio, based on M1/M0, not checking adduct form.
 
    Usage
    -----
    full_list_empCpds  = json.load(open(json_empcpd))
    isopair_empCpds = filter_khipus(full_list_empCpds)

    '''
    isopair_empCpds = []
    for epd in list_empCpds:
        M0, M1 = get_M0(epd['MS1_pseudo_Spectra']), get_M1(epd['MS1_pseudo_Spectra'])
        if M0 and M1:
            if float(M1['representative_intensity'])/(1 + float(M0['representative_intensity'])) < natural_ratio_limit:
                isopair_empCpds.append( epd['interim_id'] )

    return isopair_empCpds
    

def filter_khipus_by_samples(kdict, unlabelled=[1, 2, 3], labeled=[0, 4, 5], natural_ratio_limit=0.5):
    '''
    Enumerate khipus with good natural ratio and increased isotope labeling ratio.
    
    kdict : khipu_dict as input
    
    returns 
    list of khipus with good natural 13C ratio, list of khipus with increased 13C ratio.

    '''
    khipus_good_natural_ratios, khipus_increased_ratios = [], []
    for k,v in kdict.items():
        # interim_id = v['interim_id']
        M0, M1, Mx = get_M0(v['MS1_pseudo_Spectra']), get_M1(v['MS1_pseudo_Spectra']
                                            ), get_highest_13C(v['MS1_pseudo_Spectra'])
        if M0 and M1:
            try:
                unlabelled_13C = np.array(M1['intensities'])[unlabelled].mean()
                unlabelled_12C = np.array(M0['intensities'])[unlabelled].mean()
                ratio_natural = unlabelled_13C/(unlabelled_12C+0.1)
                if ratio_natural < natural_ratio_limit:
                    khipus_good_natural_ratios.append( (k, ratio_natural) )
                    if Mx:
                        labelled_13C = np.array(Mx['intensities'])[labeled].mean()
                        labelled_12C = np.array(M0['intensities'])[labeled].mean()
                        ratio_labeled = labelled_13C/(labelled_12C+0.1)
                        if ratio_labeled > ratio_natural:
                            khipus_increased_ratios.append( (k, ratio_labeled) )
            except IndexError:
                pass
    print(len(khipus_good_natural_ratios), len(khipus_increased_ratios))
    return khipus_good_natural_ratios, khipus_increased_ratios
    
    
    
def get_isopairs_good_khipus(list_empCpds, natural_ratio_limit=0.5):
    '''
    returns 
    Two lists of khipus, isopair_empCpds_ids (IDs only), good_khipus (full dict), and number of isopair_mtracks.
    isopair_empCpds = with good natural 13C ratio, based on M1/M0, not checking adduct form.
    good_khipus = isopair_empCpds and M0 being a good feature.
    
    Some inline MS/MS expts cause split MS1 peaks. Thus isopair_mtracks are more indicative of the data coverage.
    
    Usage
    -----
    full_list_empCpds  = json.load(open(json_empcpd))
    isopair_empCpds, num_isopair_mtracks, good_khipus = get_isopairs_good_khipus(full_list_empCpds)

    Use filter_khipus if not considering is_good_peak.
    '''
    isopair_empCpds_ids, isopair_mtracks, good_khipus = [], [], []
    for epd in list_empCpds:
        if len(epd['MS1_pseudo_Spectra']) > 1:
            try:
                M0, M1 = get_M0(epd['MS1_pseudo_Spectra']), get_M1(epd['MS1_pseudo_Spectra'])
                feature_M0 = [f for f in epd['MS1_pseudo_Spectra'] if f['isotope']=='M0'][0]
                if M0 and M1:
                    if float(M1['representative_intensity'])/(1 + float(M0['representative_intensity'])) < natural_ratio_limit:
                        isopair_empCpds_ids.append( epd['interim_id'] )
                        if feature_M0['is_good_peak']:
                            good_khipus.append( epd )
                            isopair_mtracks.append( epd["MS1_pseudo_Spectra"][0]['parent_masstrack_id'] )
            except KeyError:
                print(epd['interim_id'])
                
    return isopair_empCpds_ids, len(set(isopair_mtracks)), good_khipus
    
    
def count_singletons(list_empCpds):
    return len([epd for epd in list_empCpds if len(epd['MS1_pseudo_Spectra'])==1])
    
    
def get_isopair_features(full_list_empCpds, isopair_empCpds):
    '''
    Not clean, including more features than desired..
    '''
    isopair_features = []
    for epd in full_list_empCpds:
        if epd['interim_id'] in isopair_empCpds:
            isopair_features += epd['MS1_pseudo_Spectra']

    return isopair_features
    