'''
Class structures to connect metabolite/compound databases to empirical compounds;
and to connect experimental data to empirical compounds.

'''
import json
from mass2chem.epdsConstructor import epdsConstructor
from .search import *
from .ions import generate_ion_signature

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
                    {'feature_id': 'FT1705', 'mz': 269.0878, 'rtime': 99.90, 'charged_formula': '', 'ion_relation': 'M+H[1+]'},
                    {'feature_id': 'FT1876', 'mz': 291.0697, 'rtime': 99.53, 'charged_formula': '', 'ion_relation': 'M+Na[1+]'},
                    {'feature_id': 'FT1721', 'mz': 270.0912, 'rtime': 99.91, 'charged_formula': '', 'ion_relation': 'M(C13)+H[1+]'},
                    {'feature_id': 'FT1993', 'mz': 307.0436, 'rtime': 99.79, 'charged_formula': '', 'ion_relation': 'M+K[1+]'},
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
                   {"feature_id": peak['id'], "mz": peak['mz'], "rtime": peak['rtime'], "charged_formula": "",  "ion_relation": peak['ion_relation'],}
                )
        return {'interim_id': self.interim_id, 
                'neutral_formula_mass': self.neutral_formula_mass,
                'neutral_formula': self.neutral_formula,
                'Database_referred': self.Database_referred,
                'identity': self.write_identity(),
                'MS1_pseudo_Spectra': features,
                'MS2_Spectra': self.MS2_Spectra,
                }


class knownCompoundDatabase:
    '''
    Indexed data store for known compounds. 
    One can search by mass or mass tree (patterns of isotopes/adducts, in the form of empCpd).
    centurion_mass_tree is an indexed dictionary to look up an empCpd by 100th decimal.
    An empCpd is an empirical compound, a set of features drived from the same mass (see README), 
    which of include isomers.
    '''
    def __init__(self):
        '''

        '''
        self.mass_indexed_compounds = {}
        self.emp_cpds_trees = { 'pos': {}, 'neg': {} }           # separately for positive and negative ion modes
        
    def mass_index_list_compounds(self, list_compounds):
        '''
        Using dictionaries to index data, then compile into list_emp_cpds.
        Input
        -----
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

    def build_emp_cpds_index(self, primary_only=True):
        '''
        For each emp_cpd, generate ion signatures using C12/C13 and common_adducts adducts (pos or neg ion mode).
        Then use build_centurion_tree function from .search to build index.
        The most common ions are adequate for initial search and annoation. 
        One may choose to extend the search for other ions after the emp_cpd is matched.
        Format example -
        pos_peak_list: [{'parent_epd_id': 1670, 'mz': 133.0970237, 'ion_relation': 'M[1+]'}, ...]

        '''
        pos_peak_list, neg_peak_list = [], []
        pos_epd_ions, neg_epd_ions = {}, {}     # dictionaries using same keys as self.mass_indexed_compounds
        for k, v in self.mass_indexed_compounds.items():
            __LL = generate_ion_signature(v['neutral_formula_mass'], v['neutral_formula'], mode='pos', primary_only=True)
            # signature format e.g. [[304.203251, 'M[1+]', 'C19H28O3'], ...]
            pos_epd_ions[k] = __LL
            for ion in __LL:
                pos_peak_list.append( {'mz': ion[0], 'parent_epd_id': k, 'ion_relation': ion[1],} )
            
            # do neg ions now
            __LL = generate_ion_signature(v['neutral_formula_mass'], v['neutral_formula'], mode='neg', primary_only=True)
            neg_epd_ions[k] = __LL
            for ion in __LL:
                neg_peak_list.append( {'mz': ion[0], 'parent_epd_id': k, 'ion_relation': ion[1],} )

        # map peaks -> empCpd
        self.emp_cpds_trees['pos'] = build_centurion_tree(pos_peak_list)
        self.emp_cpds_trees['neg'] = build_centurion_tree(neg_peak_list)

    def search_mz_single(self, query_mz, mode='pos', mz_tolerance_ppm=5):
        '''
        return list of matched empCpds
        '''
        return find_all_matches_centurion_indexed_list(query_mz, self.emp_cpds_trees[mode], mz_tolerance_ppm)

    def search_mz_batch(self, query_mz_list, mode='pos', mz_tolerance_ppm=5):
        results = []
        for query_mz in query_mz_list:
            results.append(
                self.search_mz_single(query_mz, mode, mz_tolerance_ppm)
            )
        return results



    def search_emp_cpd_single(self, emp_cpd, mode='pos', mz_tolerance_ppm=5):
        pass



    def search_emp_cpd_batch(self, list_emp_cpd, mode='pos', mz_tolerance_ppm=5):
        results = []
        for emp_cpd in list_emp_cpd:
            results.append(
                self.search_emp_cpd_single(emp_cpd, self.emp_cpds_trees, mz_tolerance_ppm)
            )
        return results


    def export(self):
        pass



class ExperimentalEcpdDatabase:
    '''
    Build a boutique data store for user's experimental data, input being list of peaks or empCpds (i.e. emp_cpd).
    Allow search of Compounds, which support targeted compound search in a dataset.
    '''
    def __init__(self, mode='pos'):
        '''
        mode: ionizaation mode, 'pos' or 'neg'.
        Take input list of peaks. Peaks here are usually features.
        '''
        self.mode = mode
        self.list_peaks = []
        self.indexed_peaks = {}
        self.dict_empCpds = {}
        self.indexed_empCpds = {}

    def build_from_list_peaks(self, list_peaks):
        self.list_peaks = list_peaks
        ECCON = epdsConstructor(list_peaks, mode=self.mode)
        list_empCpds = ECCON.peaks_to_epds()
        self.dict_empCpds = self.index_reformat_epds(list_empCpds, list_peaks)
        self.index_empCpds()

    def build_from_list_empCpds(self, list_empCpds):
        for E in list_empCpds:
            self.dict_empCpds[E['interim_id']] = E
        self.index_empCpds()

    def index_empCpds(self):
        '''
        Build indices for self.list_peaks and self.empCpds.
        '''
        self.indexed_peaks  = build_centurion_tree(self.list_peaks)
        __PL = []
        for _, epd in self.dict_empCpds.items():
            peaks = epd['MS1_pseudo_Spectra']
            for P in peaks:
                P['parent_epd_id'] = epd['interim_id']
                __PL.append( P )
        self.indexed_empCpds = build_centurion_tree(__PL)

    def index_reformat_epds(self, list_empCpds, FeatureList):
        fDict = {}
        for F in FeatureList:
            fDict[F['id_number']] = F
        new = {}
        for E in list_empCpds:
            features = []
            for peak in E['list_peaks']:
                features.append(
                    {'feature_id': peak[0], 
                    'mz': fDict[peak[0]]['mz'], 
                    'rtime': fDict[peak[0]]['apex'], 
                    'charged_formula': '', 
                    'ion_relation': peak[1]}
                )
            new[E['id']] = {
                'interim_id': E['id'], 
                'neutral_formula_mass': None,
                'neutral_formula': None,
                'Database_referred': [],
                'identity': [],
                'MS1_pseudo_Spectra': features,
                'MS2_Spectra': [],
                }

        return new

    def search_peaks_mz_single(self, query_mz, mz_tolerance_ppm=5):
        '''
        return list of matched peaks
        '''
        return find_all_matches_centurion_indexed_list(query_mz, self.indexed_peaks, mz_tolerance_ppm)

    def earch_peaks_mz_batch(self, query_mz_list, mz_tolerance_ppm=5):
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



    def search_cpd_single(self, cpd, mz_tolerance_ppm=5):
        pass



    def search_cpd_batch(self, list_cpd, mz_tolerance_ppm=5):
        results = []
        for cpd in list_cpd:
            results.append(
                self.search_emp_cpd_single(cpd, self.emp_cpds_trees, mz_tolerance_ppm)
            )
        return results



    def export(self):
        '''Will add DB match too
        
        '''
        
        
        with open(outfile, 'w', encoding='utf-8') as f:
            json.dump(list_empCpds, f, ensure_ascii=False, indent=2)

        print("\nEmpirical compound annotaion (%d) was written to %s." %(len(list_empCpds), outfile))

