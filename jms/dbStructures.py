'''
Class structures to connect metabolite/compound databases to empirical compounds;
and to connect experimental data to empirical compounds.

Example compound in KCD.mass_indexed_compounds:
    ('C8H13N2O2_169.097154',
    {'interim_id': 'C8H13N2O2_169.097154',
    'neutral_formula': 'C8H13N2O2',
    'neutral_formula_mass': 169.097154087,
    'compounds': [{'primary_id': 'HMDB0062696',
        'primary_db': 'HMDB',
        'name': 'Pyridoxaminium(1+)',
        'neutral_formula': 'C8H13N2O2',
        'neutral_formula_mass': 169.097154087,
        'SMILES': 'CC1=C(O)C(C[NH3+])=C(CO)C=N1',
        'inchikey': 'NHZMQXZHNVQTQA-UHFFFAOYSA-O',
        'other_ids': {'PubChem': '25245492', 'KEGG': '', 'ChEBI': '57761'}}]})
Example search result:
        peak_result_dict can have one empCpd matched to multiple DB entries, e.g.
            {'peak': {'parent_masstrack_id': 240,
            'mz': 119.12584100277608,
            'id_number': 'F201'},
            'interim_id': 20,
            'epd_ion_relation': '13C/12C',
            'list_matches': [('C6H15NO_117.115364', 'M+H[1+]', 2),
            ('C6H13N_99.104799', 'M+H2O+H[1+]', 1)]}

'''
import json
from operator import itemgetter
import numpy as np
from mass2chem.epdsConstructor import epdsConstructor
from .search import *
from .formula import *
from .ions import compute_adducts_formulae, generate_ion_signature
from .data.list_formula_mass import list_formula_mass


def read_table_to_peaks(infile, delimiter='\t'):
    '''Merely provided as a template. Modify according to your own file format.
    '''
    list_peaks = []
    w = open(infile).readlines()
    for line in w[1:]:
        a = line.rstrip().split(delimiter)
        list_peaks.append(
            {'id_number': a[13], 'mz': float(a[2]), 
            'apex': float(a[3]), 'height': float(a[5]), 
            'cSelectivity': float(a[10]), 'goodness_fitting': float(a[11]), 'snr': float(a[12]), }
        )

    print(len(list_peaks))
    return list_peaks

def annotate_peaks_against_kcds(list_peaks, list_compounds, 
                                export_file_name_prefix='jms_annotated_',
                                mode='pos',  mz_tolerance_ppm=5):
    '''
    Wrapper function, to generate three annotation files for input list_peaks.
    list_compounds is known compound database, e.g.
    list_compounds = json.load(open('jms/data/compounds/list_compounds_HMDB4.json'))
    '''
    KCD = knownCompoundDatabase()
    KCD.mass_index_list_compounds(list_compounds)
    KCD.build_emp_cpds_index()
    KCD.export_mass_indexed_compounds(export_file_name_prefix+"KCD_mass_indexed_compounds.json")
    EED = ExperimentalEcpdDatabase()
    EED.build_from_list_peaks(list_peaks)
    search_result = EED.annotate_against_KCD(KCD)
    EED.export_annotations(search_result, KCD, export_file_name_prefix)


#----------------------------------------------------------------------------------------
class empiricalCompound:
    '''
    Simple class, following exampe in metDataModel.
    A template is:
            {
            "neutral_formula_mass": 268.08077, 
            "neutral_formula": C10H12N4O5,
            "alternative_formulas": [],
            "interim_id": C10H12N4O5_268.08077,
            "identity": [
                    {'compounds': ['HMDB0000195'], 'names': ['Inosine'], 'score': 0.6, 'probability': null},
                    {'compounds': ['HMDB0000195', 'HMDB0000481'], 'names': ['Inosine', 'Allopurinol riboside'], 'score': 0.1, 'probability': null},
                    {'compounds': ['HMDB0000481'], 'names': ['Allopurinol riboside'], 'score': 0.1, 'probability': null},
                    {'compounds': ['HMDB0003040''], 'names': ['Arabinosylhypoxanthine'], 'score': 0.05, 'probability': null},
                    ],
            "MS1_pseudo_Spectra": [
                    {'id_number': 'FT1705', 'mz': 269.0878, 'rtime': 99.90, 'charged_formula': '', 'ion_relation': 'M+H[1+]'},
                    {'id_number': 'FT1876', 'mz': 291.0697, 'rtime': 99.53, 'charged_formula': '', 'ion_relation': 'M+Na[1+]'},
                    {'id_number': 'FT1721', 'mz': 270.0912, 'rtime': 99.91, 'charged_formula': '', 'ion_relation': 'M(C13)+H[1+]'},
                    {'id_number': 'FT1993', 'mz': 307.0436, 'rtime': 99.79, 'charged_formula': '', 'ion_relation': 'M+K[1+]'},
                    ],
            "MS2_Spectra": ['AZ0000711', 'AZ0002101'],
            "Database_referred": ["Azimuth", "HMDB", "MONA"],
            }
    '''
    def __init__(self):
        self.id = id                            # e.g. 'E00001234'
        self.interim_id = ''
        self.neutral_formula_mass = None
        self.neutral_formula = ''
        self.charged_formula = ''
        self.charge = None
        self.Database_referred = []
        self.MS1_pseudo_Spectra = []            # list of features that belong to this empCpd
        self.MS2_Spectra = []                   # MS2 identifiers can be universal (e.g. hashed ids)
        self.identity = self.annotation = []    # see desired serialize() output; also in README

    def read_json_model(self, jmodel):
        self.interim_id = jmodel['interim_id']
        self.neutral_formula_mass = jmodel['neutral_formula_mass']
        self.neutral_formula = jmodel['neutral_formula']
        self.Database_referred = jmodel['Database_referred']
        self.identity = jmodel['identity']
        self.MS1_pseudo_Spectra  = jmodel['MS1_pseudo_Spectra']
        self.MS2_Spectra = jmodel['MS2_Spectra']    

    def serialize(self):
        features = []
        for peak in self.MS1_pseudo_Spectra:
                features.append(        # this is given as example; one may need to modify the mapping variable names
                   {"id_number": peak['id'], "mz": peak['mz'], "rtime": peak['rtime'], "charged_formula": "",  "ion_relation": peak['ion_relation'],}
                )
        return {'interim_id': self.interim_id, 
                'neutral_formula_mass': self.neutral_formula_mass,
                'neutral_formula': self.neutral_formula,
                'Database_referred': self.Database_referred,
                'identity': self.write_identity(),
                'MS1_pseudo_Spectra': features,
                'MS2_Spectra': self.MS2_Spectra,
                }


#----------------------------------------------------------------------------------------
class knownCompoundDatabase:
    '''
    Indexed data store for known compounds. 
    One can search by mass or mass tree (patterns of isotopes/adducts, in the form of empCpd).
    centurion_mass_tree is an indexed dictionary to group ions by 100th decimal.
    An empCpd is an empirical compound, a set of features drived from the same mass (see README), 
    which of include isomers.
    '''
    def __init__(self):
        '''
        Two main data containers, mass_indexed_compounds and emp_cpds_trees.
        The latter is indexed for searches, separately for positive and negative ion modes.
        '''
        self.mass_indexed_compounds = {}
        self.emp_cpds_trees = { 'pos': {}, 'neg': {} }
        
    def mass_index_list_compounds(self, list_compounds):
        '''
        Using dictionaries to index data, then compile into list_emp_cpds.
        Input
        =====
        list of compounds:  Example of a compound
        {
            'primary_id': HMDB0000195,
            'primary_db': 'HMDB',
            'name': 'Inosine',
            "neutral_formula": C10H12N4O5,
            "neutral_formula_mass": 268.08077, 
            'SMILES': 'OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=NC2=C(O)N=CN=C12', 
            'inchikey': 'UGQMRVRMYYASKQ-KQYNXXCUSA-N',
            'LogP': -2.10,
            'other_ids': {'PubChem': '6021',
                        'KEGG': 'C00294',
                        'ChEBI': '17596',
                        'MetaNetX': 'MNXM1103335',
                        },
        }
        See example of how to generate Compounds:
        https://github.com/shuzhao-li/JMS/wiki/Scripting-to-create-a-dictionary-from-HMDB-data

        Update
        ======
        self.mass_indexed_compounds so that isomers are grouped under same empCpd.
        '''
        _db = {}
        for cpd in list_compounds:
            k = cpd['neutral_formula']+ '_' + str(round(float(cpd['neutral_formula_mass']),6))  # ensuring unique formula and mass
            if k in _db:  
                _db[k].append( cpd )  
            else:  
                _db[k] = [cpd]  

        for k,v in _db.items():
            self.mass_indexed_compounds[k] = {
                    "interim_id": k,
                    "neutral_formula": v[0]['neutral_formula'],
                    "neutral_formula_mass": float(v[0]['neutral_formula_mass']),
                    "compounds": v,
                }

    def build_emp_cpds_index(self, primary_only=True, include_C13=False):
        '''
        For each emp_cpd, generate ion signatures common_adducts adducts (pos or neg ion mode).
        Then use build_centurion_tree function from .search to build index.

        primary_only: only considering the most common ions 
        (https://github.com/shuzhao-li/mass2chem/blob/master/mass2chem/formula.py),
        which are adequate for initial search and annoation. 
        One may choose to extend the search for other ions after the emp_cpd is matched.

        include_C13: considering 13C isotopes for the adducts.
        This is not desired in default application of matching to experimental data via empCpds.
        The isotopes and adducts should be organized into empCpds in experimental data,
        then the 12C anchor ion has to be present in search result in oder to be considered a true match.
        In that scenario, adding 13C in DB records is unecessary.
        One should include 13C when generating a database for single ion searches.

        Format example -
        pos_peak_list: [{'parent_epd_id': 1670, 'mz': 133.0970237, 'ion_relation': 'M[1+]'}, ...]

        '''
        __ion_generator__ = compute_adducts_formulae
        if include_C13:
            __ion_generator__ = generate_ion_signature

        pos_peak_list, neg_peak_list = [], []
        pos_epd_ions, neg_epd_ions = {}, {}     # dictionaries using same keys as self.mass_indexed_compounds
        for k, v in self.mass_indexed_compounds.items():
            __LL = __ion_generator__(v['neutral_formula_mass'], v['neutral_formula'], mode='pos', primary_only=True)
            # signature format e.g. [[304.203251, 'M[1+]', 'C19H28O3'], ...]
            pos_epd_ions[k] = __LL
            for ion in __LL:
                pos_peak_list.append( {'mz': ion[0], 'parent_epd_id': k, 'ion_relation': ion[1],} )
            
            # do neg ions now
            __LL = __ion_generator__(v['neutral_formula_mass'], v['neutral_formula'], mode='neg', primary_only=True)
            neg_epd_ions[k] = __LL
            for ion in __LL:
                neg_peak_list.append( {'mz': ion[0], 'parent_epd_id': k, 'ion_relation': ion[1],} )

        # map peaks -> empCpd
        self.emp_cpds_trees['pos'] = build_centurion_tree(pos_peak_list)
        self.emp_cpds_trees['neg'] = build_centurion_tree(neg_peak_list)

    def search_mz_single(self, query_mz, mode='pos', mz_tolerance_ppm=5):
        '''
        return list of matched empCpds, e.g.
            [{'mz': 130.017306555, 'parent_epd_id': 'C4H3FN2O2_130.017856', 'ion_relation': 'M[1+]'}]
        '''
        return find_all_matches_centurion_indexed_list(query_mz, self.emp_cpds_trees[mode], mz_tolerance_ppm)

    def search_mz_batch(self, query_mz_list, mode='pos', mz_tolerance_ppm=5):
        results = []
        for query_mz in query_mz_list:
            results.append(
                self.search_mz_single(query_mz, mode, mz_tolerance_ppm)
            )
        return results

    def evaluate_mass_accuracy_ratio(self, query_mz_list, mode='pos', mz_tolerance_ppm=10):
        '''
        Find the best matches of query_mz_list in known compound database, within mz_tolerance_ppm,
        then evaluate systematic shift of experimental m/z values by theoretical compound values.
        This can be used to report mass accuracy, and to help mass calibration.
        '''
        results = []        # delta m/z, [(expt - theoretical)/theoretical, ...]
        for query_mz in query_mz_list:
            found = find_best_match_centurion_indexed_list(query_mz, self.emp_cpds_trees[mode], mz_tolerance_ppm)
            if found:
                results.append( query_mz/found['mz'] - 1 )

        if results:
            ratio = np.mean(results)
            print("Mass accuracy was estimated on %d matched values as %2.1f ppm." %(len(results), ratio*1000000))
            return ratio
        else:
            print("No m/z match found in database to estimate mass accuracy.")
            return None

    def search_emp_cpd_single(self, emp_cpd, mode='pos', mz_tolerance_ppm=5):
        '''
        emp_cpd format has to follow (loosely) specificiations in metDataModel, requiring anchor ion, e.g.
        {interim_id': 12,
        'neutral_formula_mass': None,
        'neutral_formula': None,
        'MS1_pseudo_Spectra': [{'id_number': 'F94',
                                'mz': 98.97531935309317,
                                'rtime': 700.0,
                                'charged_formula': '',
                                'ion_relation': 'anchor',
                                'parent_epd_id': 12}, ], ...}
        Note that neutral_formula is not assigned from experimental measurements, which only indicate m/z distances
        between coeluting ions. 
        The emp_cpds_trees (from mass_indexed_compounds) only considers primary ions here, and
        should NOT include isotopes, which are considered in the experimental data.
        The ion_relations are expected to be ['M[1+]', 'M+H[1+]', 'M+Na[1+]', 'M+H2O+H[1+]'] for pos,
        and ['M[-]', 'M-H[-]', 'M-H2O-H[-]', 'M+Cl[-]'] for neg.

        Return
        ======
        list of matched db empCpd as [(empCpd_id, ion_relation, score), ...], sorted by score.
        While most search will return a single match, undeterministic situations exist,
        e.g. C6H13N and C6H14N could match to M+ and M+H+, respectively.
        The resolution of these will be left to downstream methods (bayesian etc).
        The score is number of matched ions.
        '''
        results = []
        query_mzs = [x['mz'] for x in emp_cpd['MS1_pseudo_Spectra']]
        # look for anchor ion first; default first ion in 'MS1_pseudo_Spectra'
        matches = self.search_mz_single(query_mzs[0], mode, mz_tolerance_ppm)
        # format - [{'mz': 130.017306555, 'parent_epd_id': 'C4H3FN2O2_130.017856', 'ion_relation': 'M[1+]'}, ...]
        for _M in matches:
            R = self.mass_indexed_compounds[_M['parent_epd_id']]
            # compute_adducts_formulae considers both isotopes and adducts
            db_peaks = compute_adducts_formulae(R['neutral_formula_mass'], R['neutral_formula'], mode)
            results.append((_M['parent_epd_id'],
                            _M['ion_relation'],     # ion_relation btw anchor and DB record, not ions in empCpd
                            score_emp_cpd_matches(query_mzs, [x[0] for x in db_peaks], mz_tolerance_ppm), 
                    ))
        return sorted(results, key=itemgetter(2), reverse=True)
        
    def search_emp_cpd_batch(self, list_emp_cpd, mode='pos', mz_tolerance_ppm=5):
        results = []
        for emp_cpd in list_emp_cpd:
            results.append(
                self.search_emp_cpd_single(emp_cpd, mode, mz_tolerance_ppm)
            )
        return results

    def short_report_emp_cpd(self, interim_id):
        return ";".join([ C['primary_id']+'$'+C['name'] 
                        for C in self.mass_indexed_compounds[interim_id]['compounds'] ])

    def export_mass_indexed_compounds(self, outfile="KCD_mass_indexed_compounds.json"):
        with open(outfile, 'w', encoding='utf-8') as O:
            json.dump(self.mass_indexed_compounds, O, ensure_ascii=False, indent=2)

    def export_search_emp_cpd_batch(self, batch_result):
        '''
        '''
        pass


#----------------------------------------------------------------------------------------
class ExperimentalEcpdDatabase:
    '''
    Build a boutique data store for user's experimental data, input being list of peaks or empCpds (i.e. emp_cpd).
    Allow search of Compounds, which support targeted compound search in a dataset.

    The construction of empirical compounds (empCpds) is based on 
    1) initial search of ion patterns (C13, +H, +Na, -H, +Na-2H +Cl-) based on mass2chem.epdsConstructor
    2) match to databased generated through knownCompoundDatabase
    3) augment formula match via data.list_formula_mass (HMDB4+PubChemLite)
    4) formula based grid search to extend all empCpds

    self.dict_empCpds is updated in situ.
    '''
    def __init__(self, mode='pos'):
        '''
        mode: ionizaation mode, 'pos' or 'neg'.
        Take input list of peaks. Peaks here are usually features.
        '''
        self.mode = mode
        self.list_peaks = []
        self.indexed_peaks = {}
        self.dict_peaks = {}
        self.dict_empCpds = {}
        self.indexed_empCpds = {}
        self.peak_to_empCpd = {}
        self.peak_to_empCpd_ion_relation = {}

    def build_from_list_peaks(self, list_peaks):
        '''
        list of peaks, e.g. [ {'id_number': 555,        # change to 'id_number' throughout
                                'mz': 133.0970, 
                                'apex': 654, 
                                'height': 14388.0, ...
                                }, ...]
        '''
        self.list_peaks = list_peaks
        ECCON = epdsConstructor(list_peaks, mode=self.mode)
        self.dict_empCpds = ECCON.peaks_to_epdDict(
                seed_search_patterns = ECCON.seed_search_patterns, 
                ext_search_patterns = ECCON.ext_search_patterns,
                mz_tolerance_ppm=5, 
                coelution_function='overlap',
                check_isotope_ratio = True
        ) 
        self.index_empCpds()

    def build_from_list_empCpds(self, list_empCpds):
        for E in list_empCpds:
            self.dict_empCpds[E['interim_id']] = E
        self.index_empCpds()

    def index_empCpds(self):
        '''
        Build indices for self.list_peaks and self.empCpds.
        Minor possibility that peak:empCpd is not N:1.
        round 1 index of self.dict_empCpds; round 2 to be done after annotation
        '''
        for P in self.list_peaks:
            self.dict_peaks[P['id_number']] = P

        self.indexed_peaks  = build_centurion_tree(self.list_peaks)
        self.formula_tree = get_formula_ions_tree(list_formula_mass, mode=self.mode)

        __PL = []
        for interim_id, epd in self.dict_empCpds.items():
            peaks = epd['MS1_pseudo_Spectra']
            for P in peaks:
                self.peak_to_empCpd[P['id_number']] = interim_id
                self.peak_to_empCpd_ion_relation[P['id_number']] = P['ion_relation']
                P['parent_epd_id'] = epd['interim_id'] 
                __PL.append( P )

        self.indexed_empCpds = build_centurion_tree(__PL)

    # standalone search functions
    def search_peaks_mz_single(self, query_mz, mz_tolerance_ppm=5):
        '''
        return list of matched peaks
        '''
        return find_all_matches_centurion_indexed_list(query_mz, self.indexed_peaks, mz_tolerance_ppm)

    def search_peaks_mz_batch(self, query_mz_list, mz_tolerance_ppm=5):
        results = []
        for query_mz in query_mz_list:
            results.append(
                self.search_peaks_mz_single(query_mz, mz_tolerance_ppm)
            )
        return results

    def search_empCpds_mz_single(self, query_mz, mz_tolerance_ppm=5):
        '''
        return list of matched empCpds
        '''
        return find_all_matches_centurion_indexed_list(query_mz, self.indexed_empCpds, mz_tolerance_ppm)

    def search_empCpds_mz_batch(self, query_mz_list, mz_tolerance_ppm=5):
        results = []
        for query_mz in query_mz_list:
            results.append(
                self.search_empCpds_mz_single(query_mz, mz_tolerance_ppm)
            )
        return results

    def search_peaks_compound_single(self, compound, mz_tolerance_ppm=5):
        '''
        Search a compound against this experimental dataset, and return matched peaks.
        Compound format:
        { 'primary_id': HMDB0000195, 'name': 'Inosine',
            "neutral_formula": C10H12N4O5,
            "neutral_formula_mass": 268.08077, ..}
        '''
        expected_peaks = compute_adducts_formulae(
            compound['neutral_formula_mass'], compound['neutral_formula'], self.mode, primary_only=False)
        # format [[304.203251, 'M[1+]', 'C19H28O3'], ...]
        return list(zip( [x[1] for x in expected_peaks],
                    self.search_peaks_mz_batch([x[0] for x in expected_peaks], mz_tolerance_ppm)
                    ))

    def search_peaks_compound_batch(self, list_cpd, mz_tolerance_ppm=5):
        results = []
        for cpd in list_cpd:
            results.append(
                self.search_peaks_compound_single(cpd, cpd, mz_tolerance_ppm)
            )
        return results

    def annotate_empCpds_against_KCD(self, KCD):
        '''
        Get all empCpd matches between this experimental dataset and a known compound database.
        KCD: knownCompoundDatabase instance.
        >>> KCD.search_emp_cpd_single( EED.dict_empCpds[15] )
        [('C6H14N_100.112624', 'M[1+]', 2), ('C6H13N_99.104799', 'M+H[1+]', 2)]
        This search empCpds first, in which only anchor ion is searched. 
        Because if anchor ion is not found in the database, other adducts are not possible.

        return dictionary e.g.
                {(101, []),
                (102,
                [('C9H20NO2_174.149404', 'M[1+]', 2),
                ('C9H17NO_155.131014', 'M+H2O+H[1+]', 1)]),
                (103,
                [('C6H15N4O2_175.119501', 'M[1+]', 2),
                ('C6H14N4O2_174.111676', 'M+H[1+]', 2)]),
                (104, []), ...}
        '''
        resultDict = {}
        for epd in self.dict_empCpds.values():
            resultDict[epd['interim_id']] = KCD.search_emp_cpd_single(epd, self.mode)
        return resultDict

    def search_mz_for_formula(self, mz):
        '''
        return best matched formula or None, 
        format as {'mz': ion[0], 'neutral_formula': formula, 
                               'ion_relation': ion[1],
                               'neutral_formula_mass': mass}
        '''
        return search_mz_formula_tree(mz, self.formula_tree, limit_ppm=5)

    # Annotation functions
    def extend_empCpd_annotation(self, KCD):
        self.empCpds_formula_search(KCD)
        self.__extend_empCpds__()

    def empCpds_formula_search(self, KCD):
        '''
        Update self.dict_empCpds for formulae first by KCD search then .data.formula_tree.
        KCD: knownCompoundDatabase instance.
        '''
        epd_search_result_dict = self.annotate_empCpds_against_KCD(KCD)
        for interim_id, V in epd_search_result_dict.items():
            if V:
                # update empCpd with formula matches
                self.dict_empCpds[interim_id]['list_matches'] = V
                self.dict_empCpds[interim_id]['neutral_formula'] = \
                    KCD.mass_indexed_compounds[V[0][0]]['neutral_formula']
                self.dict_empCpds[interim_id]['neutral_formula_mass'] = \
                    KCD.mass_indexed_compounds[V[0][0]]['neutral_formula_mass']
            else:
                # first peak should be anchor
                anchor_mz = self.dict_empCpds[interim_id]['MS1_pseudo_Spectra'][0]['mz']
                F = self.search_mz_for_formula(anchor_mz)
                if F:
                    self.dict_empCpds[interim_id]['neutral_formula'] = F['neutral_formula']
                    self.dict_empCpds[interim_id]['neutral_formula_mass'] = F['neutral_formula_mass']

    def __extend_empCpds__(self, mz_tolerance_ppm=5):
        '''
        Extend empCpds by 
        additional isotopes/adducts, based on mass2chem.formula.compute_adducts_formulae
        '''
        peakList = [P for P in self.list_peaks if P['id_number'] not in self.peak_to_empCpd]
        peakTree = build_centurion_tree(peakList)
        for k,E in self.dict_empCpds.items():
            if E['neutral_formula']:
                anchor = E['MS1_pseudo_Spectra'][0]
                expected_ions = compute_adducts_formulae(E['neutral_formula_mass'], 
                                                E['neutral_formula'], self.mode, primary_only=False)
                                                # [(58.53894096677, 'M+2H[2+]', result_formula), ...,]
                for ion in expected_ions:
                    matches = find_all_matches_centurion_indexed_list(ion[0], peakTree, mz_tolerance_ppm)
                    for peak in matches:
                        if is_coeluted(anchor, peak, rt_tolerance=10):
                            peak['ion_relation'] = ion[1]
                            E['MS1_pseudo_Spectra'].append(peak)
                            self.peak_to_empCpd_ion_relation[peak['id_number']] = peak['ion_relation']
                            self.peak_to_empCpd[peak['id_number']] = E['interim_id']
                            # this may overwrite 1:N relationships

    def singleton_formula_search(self, KCD):
        '''
        Search singletons for formulae first by KCD search then .data.formula_tree.
        KCD: knownCompoundDatabase instance.
        '''
        found = []
        singletons = [p for p in self.dict_peaks if p not in self.peak_to_empCpd.keys()]
        for p in singletons:
            _mz = self.dict_peaks[p]['mz']
            list_matches = KCD.search_mz_single(_mz)
            # [{'mz': 130.017306555, 'parent_epd_id': 'C4H3FN2O2_130.017856', 'ion_relation': 'M[1+]'}]
            if list_matches:
                # take 1st match only here
                _epd = KCD.mass_indexed_compounds[list_matches[0]['parent_epd_id']]
                found.append((p, _epd))
            else:
                # formula search
                F = self.search_mz_for_formula(_mz)
                if F:
                    found.append((p, F))            # p is id_number

        return found

    def annotate_singletons(self, KCD):
        '''
        Search singletons for formulae first by KCD search then .data.formula_tree.
        Add new empCpds to self.dict_empCpds as new empCpds, and extend adduct search.
        '''
        formula_to_peaks = {}
        found = self.singleton_formula_search(KCD)
        for p,F in found:
            k = F['neutral_formula']
            if k in formula_to_peaks:
                formula_to_peaks[k].append((p, F))
            else:
                formula_to_peaks[k] = [(p, F)]

        found_peaks = set([x[0] for x in found])
        peakList = [P for P in self.list_peaks if P['id_number'] not in self.peak_to_empCpd and
                                                  P['id_number'] not in found_peaks]
        peakTree = build_centurion_tree(peakList)

        new_id_start = max(self.dict_empCpds.keys())
        for formula, PP in formula_to_peaks.items():
            neutral_formula_mass = PP[0][1]['neutral_formula_mass']
            P1 = self.dict_peaks[PP[0][0]]
            tmp = [P1, ]
            for jj in range(1, len(PP)):
                _P = self.dict_peaks[PP[jj][0]]
                if is_coeluted(tmp[-1], _P, rt_tolerance=10):
                    tmp.append(_P)
                else:
                    # not coeluted, new empCpd
                    new_id_start += 1
                    self.dict_empCpds[new_id_start] = {'interim_id': new_id_start,
                            'neutral_formula_mass': neutral_formula_mass, 'neutral_formula': formula,
                            'MS1_pseudo_Spectra': self.__extend_peakList__(
                                formula, neutral_formula_mass, tmp, peakTree, mz_tolerance_ppm=5),
                    }
                    tmp = [_P, ]

            new_id_start += 1
            self.dict_empCpds[new_id_start] = {'interim_id': new_id_start,
                    'neutral_formula_mass': neutral_formula_mass, 'neutral_formula': formula,
                    'MS1_pseudo_Spectra': self.__extend_peakList__(
                                formula, neutral_formula_mass, tmp, peakTree, mz_tolerance_ppm=5),
            }


    def __extend_peakList__(self, neutral_formula, neutral_formula_mass, epd_peaks, peakTree,
            mz_tolerance_ppm=5):
        '''
        Extend epd_peaks by 
        additional isotopes/adducts, based on mass2chem.formula.compute_adducts_formulae, which 
        produces [(58.53894096677, 'M+2H[2+]', result_formula), ...,]
        '''
        new = []
        anchor = epd_peaks[0]
        for ion in compute_adducts_formulae(neutral_formula_mass, neutral_formula,
                                        self.mode, primary_only=False):
            matches = find_all_matches_centurion_indexed_list(ion[0], peakTree, mz_tolerance_ppm)
            for peak in matches:
                if is_coeluted(anchor, peak, rt_tolerance=10):
                    peak['ion_relation'] = ion[1]
                    new.append(peak)
        return epd_peaks + new


    def annotate_all_against_KCD(self, KCD):
        '''
        Get matches of both empCpds and singleton peaks only in KCD.
        This function should not be used if one uses empCpds_formula_search and annotate_by_formula_grid.
        Results contain annotation via empCpd, via singleton, or empCpd features as if singleton.

        KCD: knownCompoundDatabase instance.
        return: peak_result_dict, epd_search_result_dict
        '''
        peak_result_dict = {}
        epd_search_result_dict = self.annotate_empCpds_against_KCD(KCD)
        for peak in self.list_peaks:
            ii, mz, apex = peak['id_number'], peak['mz'], peak['apex']
            matched_DB_shorts, matched_DB_records = '', ''
            interim_id = self.peak_to_empCpd.get(ii, None)  # id on expt data; ii could be 0
            epd_ion_relation = ''
            list_matches = []
            if interim_id != None:
                epd_ion_relation = self.peak_to_empCpd_ion_relation[ii]
                # Meaning there's empCpd for this feature, we look for earlier search result in search_result_dict
                list_matches = epd_search_result_dict[interim_id]
                # [('C7H8N2O_136.063663', 'M+H[1+]', 2), ('C7H9N2O_137.071488', 'M[1+]', 2), ('C7H6N2_118.053098', 'M+H2O+H[1+]', 1)]
                # this can be empty, but means the empCpd is not found in KCD.
            if not list_matches:
                # Meaning singleton feature OR unmatched empCpd feature, new search by m/z
                list_matches = KCD.search_mz_single(mz)
                # [{'mz': 130.017306555, 'parent_epd_id': 'C4H3FN2O2_130.017856', 'ion_relation': 'M[1+]'}]
                list_matches = [(x['parent_epd_id'], x['ion_relation'], None) for x in list_matches]
            if list_matches:
                peak_result_dict[ii] = {'peak': peak,
                                    'interim_id': interim_id,               # this is empCpd id
                                    'epd_ion_relation': epd_ion_relation,
                                    'list_matches': list_matches,
                                    }

        return peak_result_dict, epd_search_result_dict


    def choose_top_epd(list_KCD_matches):
        '''recommend best from matches, e.g. 
        [('C6H14N_100.112624', 'M[1+]', 2), ('C6H13N_99.104799', 'M+H[1+]', 2)]
        by 1) highest score, and 2) peck order of ions (M+H > M+Na > M+ > M+H2O+H)
        '''
        pass


    def export_empCpds(self, outfile="EED_empCpds.json"):
        with open(outfile, 'w', encoding='utf-8') as O:
            json.dump(list(self.dict_empCpds.values()), O, ensure_ascii=False, indent=2)


    def export_annotations(self, KCD, export_file_name_prefix, include_succinct=False):
        '''
        Export expt empCpds, match relationships, and a tsv table for input list_peaks, with all matched db empCpds.

        To-do:
        Will add option to include a version succinct annotation of anchor ions and recommended matches only.
        '''
        peak_result_dict, epd_search_result_dict = self.annotate_all_against_KCD(KCD)

        self.export_empCpds(export_file_name_prefix+"EED_empCpds.json")
        with open(export_file_name_prefix+"mapping.json", 'w', encoding='utf-8') as O:
            json.dump(epd_search_result_dict, O, ensure_ascii=False, indent=2)

        s = "[peak]id_number\tmz\tapex\t[EmpCpd]interim_id\t[EmpCpd]ion_relation\tmatched_DB_shorts\tmatched_DB_records\n"
        for ii, V in peak_result_dict.items():
            matched_DB_shorts, matched_DB_records = '', ''
            list_matches = V['list_matches']
            if list_matches:
                    matched_DB_shorts = ", ".join([ "(" + KCD.short_report_emp_cpd(xx[0]) + ")"  for xx in list_matches])
                    matched_DB_records = ", ".join([str(xx) for xx  in list_matches])

            s += '\t'.join([str(x) for x in [
                ii, V['peak']['mz'], V['peak']['apex'], V['interim_id'], V['epd_ion_relation'],
                matched_DB_shorts, matched_DB_records]]) + "\n"

        outfile = export_file_name_prefix + '.tsv'
        with open(outfile, 'w') as O:
            O.write(s)

        print("\nAnnotation of %d Empirical compounds was written to %s." %(len(self.dict_empCpds), outfile))


    #-----------------------------------------

    def update_annotation(self, KCD, peak_result_dict, adduct_patterns):
        '''
        Update self.dict_empCpds, and extend search by formula-based additional isotopes/adducts.
        Not updating indexing or peak to empCpd mapping.

        peak_result_dict: result from self.annotate_all_against_KCD, which
         contains annotation via empCpd, via singleton, or empCpd features as if singleton.


        '''
        
        cpd_index_to_peaks = {}
        for p,v in peak_result_dict.items():
            cpd_index = v['list_matches'][0][0]
            if v['interim_id']:
                # peak has empCpd, update empCpd with formula matches
                self.dict_empCpds[v['interim_id']]['list_matches'] = v['list_matches']
                self.dict_empCpds[v['interim_id']]['neutral_formula'] = \
                    KCD.mass_indexed_compounds[cpd_index]['neutral_formula']
                self.dict_empCpds[v['interim_id']]['neutral_formula_mass'] = \
                    KCD.mass_indexed_compounds[cpd_index]['neutral_formula_mass']
            else:
                # compile to cpd_index_to_peaks, in order to group peaks by cpd_index e.g. `C6H15NO3_149.105193`
                if cpd_index in cpd_index_to_peaks:
                    cpd_index_to_peaks[cpd_index].append(p)
                else:
                    cpd_index_to_peaks[cpd_index] = [p]

        # create new empCpds based on cpd_index_to_peaks
        new_id_start = max(self.dict_empCpds.keys())
        for cpd_index, peaks in cpd_index_to_peaks.items():
            new_id_start += 1
            self.dict_empCpds[new_id_start] = {'interim_id': new_id_start,
                    'neutral_formula_mass': KCD.mass_indexed_compounds[cpd_index]['neutral_formula_mass'],
                    'neutral_formula': KCD.mass_indexed_compounds[cpd_index]['neutral_formula'],
                    'MS1_pseudo_Spectra': [self.dict_peaks[p] for p in peaks]
            }

        # try to assign remaining peaks to empCpds
        remaining_peaks = [p for p in self.dict_peaks if p not in 
                            set(list(peak_result_dict.keys()) + list(self.peak_to_empCpd.keys()))]

        ECCON = epdsConstructor(remaining_peaks, mode=self.mode)
        list_empCpds = ECCON.extend_empCpds_by_formula_grid(

            self.dict_empCpds.values(), remaining_peaks, adduct_patterns, mz_tolerance_ppm=5,


        ) 
        # update self.dict_empCpds
        for EC in list_empCpds:
            self.dict_empCpds[EC['interim_id']] = EC
