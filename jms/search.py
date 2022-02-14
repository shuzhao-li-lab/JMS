'''

Following tools in mass2chem.

tree = build_centurion_tree(list_peaks)

find_all_matches_centurion_indexed_list(112.0950, tree, limit_ppm=5)
Out[17]:
[{'id_number': 'F241',
  'mz': 112.09503429766107,
  'apex': '267',
  'height': '3217405.650383509',
  'cSelectivity': '0.6596858638743456',
  'goodness_fitting': '0.983822971450631',
  'snr': '4'},
 {'id_number': 'F242',
  'mz': 112.09503429766107,
  'apex': '408',
  'height': '2461923.7027675146',
  'cSelectivity': '0.1015625',
  'goodness_fitting': '0.7678341743874929',
  'snr': '3'}]


Handling degenerate peaks in experimental data vs theoretical compounds 
involves two opposite problems.
We are not using the identical set of isotopes and ions yet -
they have different computational complexities, and it's not clear that we should. 

'''
isotopic_patterns = [
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

common_adducts = {
    # mass diff, modification
    # not using (intensity ratio constraint), but it can be documented or learned
    'pos': [
        (1.0078, 'H'),
        (21.9820, 'Na/H'), # Na replacing H
        (10.991, 'Na/H, double charged'),
        (18.0106, '+H2O'), 
        (18.033823, '+NH4'),
        (37.9559, '39K/H'),
        (39.9540, '41K/H'),
        (41.026549, 'Acetonitrile'),
    ],
    'neg': [
        (1.0078, 'H'),
        (22.9893, 'Na'),
        (20.97474706646, '+Na-2H'),
        (18.0106, 'H2O'), 
        (34.9689, '35Cl'),
        (36.9659, '37Cl'),
        (40.01926853323, '+ACN-H'),
        (44.998201, 'COOH'),
        (59.013295, 'CH3COO'),
    ],
}

extended_adducts = {
    'pos': [],
    'neg': [],
}


#
# -----------------------------------------------------------------------------
#

def build_centurion_tree(list_peaks):
    '''
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
                    'left_base': 648, 'right_base': 655, 'id_number': 555}, ...]
    Return a dictionary, indexing mzList by 100*mz bins.
    Because most high-resolution mass spectrometers measure well under 0.01 amu, 
    one only needs to search the corresponding 0.01 bin and two adjacent bins (to capture bordering values).
    list_mass_tracks has similar format as list_peaks.
    '''
    d = {}
    for p in list_peaks:
        cent = int(100 * p['mz'])
        if cent in d:
            d[cent].append(p)
        else:
            d[cent] = [p]
    return d


def build_peak_id_dict(list_peaks):
    d = {}
    for p in list_peaks:
        d[p['id_number']] = p
    return d


def __build_centurion_tree_mzlist(mzList):
    '''
    Return a dictionary, indexing mzList by 100*mz bins.
    Because most high-resolution mass spectrometers measure well under 0.01 amu, 
    one only needs to search the corresponding 0.01 bin and two adjacent bins (to capture bordering values).
    '''
    d = {}
    for ii in range(len(mzList)):
        cent = int(100*mzList[ii])
        if cent in d:
            d[cent].append((mzList[ii], ii))
        else:
            d[cent] = [(mzList[ii], ii)]
    return d


def find_all_matches_centurion_indexed_list(query_mz, mz_centurion_tree, limit_ppm=5):
    '''
    Return matched peaks in mz_centurion_tree.
    '''
    q = int(query_mz * 100)
    mz_tol = query_mz * limit_ppm * 0.000001
    results = []
    for ii in (q-1, q, q+1):
        L = mz_centurion_tree.get(ii, [])
        for peak in L:
            if abs(peak['mz']-query_mz) < mz_tol:
                results.append(peak)
                
    return results


def find_best_match_centurion_indexed_list(query_mz, mz_centurion_tree, limit_ppm=2):
    '''
    Return matched indices in mz_centurion_tree (based on peak list).
    '''
    q = int(query_mz * 100)
    mz_tol = query_mz * limit_ppm * 0.000001
    result = (None, 999)
    for ii in (q-1, q, q+1):
        L = mz_centurion_tree.get(ii, [])
        for peak in L:
            _d = abs(peak['mz']-query_mz)
            if _d < min(result[1], mz_tol):     # enforce mz_tol here
                result = (peak, _d)
                
    return result[0]


def is_coeluted(P1, P2, rt_tolerance=10):
    '''
    coelution is defined by overlap more than half of the smaller peak, or apexes within rt_tolerance.
    
    Example peak format: {'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0,
    'left_base': 648, 'right_base': 655, 'id_number': 555}
    rt_tolerance: tolerance of retention time is used when peak boundaries are not given.
                Can be scan numbers or seconds.
                
    return: True or False
    '''
    _coeluted = False
    try:
        len1, len2 = P1['right_base'] - P1['left_base'], P2['right_base'] - P2['left_base']
        # overlap is the max L and min R
        overlap = min(P1['right_base'], P2['right_base']) - max(P1['left_base'], P2['left_base'])
        if overlap > 0.5 * min(len1, len2):
            _coeluted = True
            
    except KeyError:
        if abs(P1['apex']-P2['apex']) <= rt_tolerance:
            _coeluted = True
        
    return _coeluted

   
def find_isotopic_signatures(list_peaks, mztree, isotopic_patterns, mz_tolerance_ppm=5, rt_tolerance_scans=5):
    '''
    See find_isotopic_pairs. This extends to all related isotopic signatures.
    Allows ambiguous matches, which are dealt with in empCpd statistics.
    Input
    =====
    isotopic_patterns = [(1.003355, '13C/12C', (0, 0.2)), ...], the third item is optional limits of abundance ratio.
    Return
    ======
    list of lists of peak numbers that match isotopic patterns, [ [], ... ]
    Example
    =======
    [ [(195, 'anchor'), (206, '13C/12C')], 
      [(182, 'anchor'), (191, '13C/12C'), (205, '18O/16O')],
      [(295, 'anchor'), (335, '13C/12C'), (368, 'M(13C),M(34S)')], ...]
    '''

    signatures = []
    for P1 in list_peaks:
        matched = [ (P1['id_number'], 'anchor'), ]          # Nature is nice to have lowest mass for the most abundant 
        for isotopic_pair in isotopic_patterns:
            (mass_difference, relation) = isotopic_pair[:2]
            tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, mztree, mz_tolerance_ppm)
            for P2 in tmp:
                if abs(P1['apex']-P2['apex']) <= rt_tolerance_scans and is_coeluted(P1, P2):
                    if len(isotopic_pair) == 2:             # not checking abundance ratio
                        matched.append( (P2['id_number'], relation) )
                    else:
                        (abundance_ratio_min, abundance_ratio_max) = isotopic_pair[2]
                        # checking abundance ratio
                        if abundance_ratio_min*P1['height'] < P2['height'] < abundance_ratio_max*P1['height']:
                            matched.append( (P2['id_number'], relation) )
        if len(matched) > 1:
            signatures.append(matched)

    return signatures


def find_adduct_signatures(list_peaks, mztree, adduct_patterns, mz_tolerance_ppm=5):
    '''
    Search adduct mass_diff in ceelution peaks, de novo. 
    Not requiring matched apex as in isotopic pairs, nor abundance ratio restriction.
    '''
    signatures = []
    for P1 in list_peaks:
        matched = [ (P1['id_number'], 'anchor'), ]
        for adduct in adduct_patterns:
            tmp = find_all_matches_centurion_indexed_list( P1['mz'] + adduct[0], mztree, mz_tolerance_ppm )
            for P2 in tmp:
                if is_coeluted(P1, P2):
                    matched.append( (P2['id_number'], adduct[1]) )

        if len(matched) > 1:
            signatures.append(matched)

    return signatures


def find_mzdiff_pairs_from_masstracks(list_mass_tracks, list_mz_diff=[1.003355, 21.9820], mz_tolerance_ppm=5):
    '''
    Find all pairs in list_mass_tracks that match a pattern in list_mz_diff, and return their id_numbers as pairs.
    This function does not use coeluction (rtime) rules. 
    Input
    =====
    list_mass_tracks: [{ 'id_number': ii,  'mz': xic[0], 'rt_scan_numbers': xic[1],  'intensity': xic[2],  }, ...]
    list_mz_diff: defaul 1.003355, 21.9820 are 13C/12C, Na/H. Dependent on ionization mode.
    
    Return
    ======
    list of pairs of mass tracks numbers.
    '''
    pairs = []
    # list_mass_tracks has similar format as list_peaks.
    mztree = build_centurion_tree(list_mass_tracks)
    for mzdiff in list_mz_diff:
        for P1 in list_mass_tracks:
            P2 = find_best_match_centurion_indexed_list(P1['mz'] + mzdiff, mztree, mz_tolerance_ppm)
            if P2:
                pairs.append((P1['id_number'], P2['id_number']))

    return pairs

def score_emp_cpd_matches(L1_peaks, L2_peaks, ppm):
    '''Return number of matched peaks
    '''
    return len(mass_paired_mapping(L1_peaks, L2_peaks, ppm)[0])


#------------------------------------------------------------------------------------------------------------------------------
# From asari.mass_functions

def mass_paired_mapping(list1, list2, std_ppm=5):
    '''
    To find unambiguous matches of m/z values between two lists.
    This sorts all m/z values first, then compare their differences in sequential neighbors.
    To be considered as an unambiguous match, the m/z values from two lists 
    should have no overlap neighbors in either direction in either list other than their own pair.
    Not necessary to have full mapping. FeatureMap is done in multiple steps. This is step 1.
    This shares some similarity to the RANSAC algorithm but prioritizes selectivity.
    For illustration, one can use one-step Gaussian model for mass shift.
    Since only mean shift is used here, and stdev is implicitly enforced in matching, no need to do model fitting.
    Input
    =====
    Two lists of m/z values, not ncessarily same length.
    std_ppm: instrument accuracy to guide value matching. Low-selectiviy values are not considered in matching.
    Return
    ======
    mapped: mapping list [(index from list1, index from list2), ...]
    ratio_deltas: mean m/z ratio shift between two lists. This is ppm*10^-6. No need to convert btw ppm here.
    Test
    ====
    list1 = [101.0596, 101.061, 101.0708, 101.0708, 101.1072, 101.1072, 101.1072, 102.0337, 102.0337, 102.0548, 102.0661, 102.0912, 102.0912, 102.1276, 102.1276, 103.0501, 103.0501, 103.0541, 103.0865, 103.0865, 103.9554, 104.0368, 104.0705, 104.0705, 104.1069, 104.1069, 104.9922, 105.0422, 105.0698, 105.0698, 105.0738, 105.1039, 105.1102, 105.9955, 106.0497, 106.065, 106.065, 106.0683, 106.0683, 106.0861, 106.0861, 106.0861, 106.1111, 106.9964, 107.0475, 107.0602, 107.0653, 107.0895, 107.9667, 108.0443, 108.0555, 108.0807, 109.0632, 109.0759]
    list2 = [101.0087, 101.035, 101.0601, 101.0601, 101.0601, 101.0601, 101.0713, 101.0714, 101.1077, 101.1077, 101.1077, 101.1077, 101.1077, 101.1158, 101.1158, 102.0286, 102.0376, 102.0468, 102.0539, 102.0554, 102.0554, 102.0554, 102.0554, 102.0666, 102.0917, 102.0917, 102.0917, 102.0918, 102.1281, 102.1281, 102.1282, 103.0394, 103.0505, 103.0507, 103.0547, 103.1233, 103.8162, 103.956, 103.956, 103.956, 104.0532, 104.0533, 104.0641, 104.0709, 104.071, 104.0831, 104.0878, 104.0895, 104.0953, 104.1073, 104.1073, 104.1074, 104.1074, 104.1182, 104.1199, 104.1265, 104.1318, 104.1354, 104.1725, 104.3998, 104.9927, 104.9927, 104.9927, 104.9927, 105.0654, 105.0703, 105.1043, 105.1133, 106.049, 106.0503, 106.0655, 106.0688, 106.0866, 106.0867, 106.0867, 106.0867, 106.114, 107.048, 107.0481, 107.0496, 107.0608, 107.0658, 108.0109, 108.0482, 108.0604, 108.0812, 108.0812, 108.9618, 109.0507, 109.0637, 109.0637, 109.0764, 109.1015]
    mass_paired_mapping(list1, list2) >>>
        ([(10, 23), (29, 65), (31, 66), (36, 70), (38, 71), (46, 81), (53, 91)],
        [4.898762180656323e-06,
        4.758718686464085e-06,
        3.805743437700149e-06,
        4.714068193732999e-06,
        4.713921530199148e-06,
        4.670025348919892e-06,
        4.583942997773922e-06])
    '''
    all = [(list1[ii], 1, ii) for ii in range(len(list1))] + [(list2[jj], 2, jj) for jj in range(len(list2))]
    # [(mz, list_origin, index_origin), ...]
    all.sort()
    NN = len(all)
    # Add a mock entry to allow loop goes through NN.
    all.append((999999, 2, None))
    mapped, ratio_deltas = [], []
    for ii in range(1, NN):
        if all[ii][1] != all[ii-1][1]:          # from two diff list_origin
            _tolerance = all[ii][0] * std_ppm * 0.000001
            _d = all[ii][0]-all[ii-1][0]
            if _d < _tolerance and all[ii+1][0]-all[ii][0] > _tolerance:
                # not allowing ii to be matched to both ii-1 and ii+1
                if all[ii][1] > all[ii-1][1]:   # always ordered as list1, list2
                    mapped.append( (all[ii-1][2], all[ii][2]) )
                    ratio_deltas.append( _d/all[ii][0] )
                else:
                    mapped.append( (all[ii][2], all[ii-1][2]) )
                    ratio_deltas.append( -_d/all[ii][0] )

    return mapped, ratio_deltas
