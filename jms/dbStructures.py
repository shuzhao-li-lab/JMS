'''
Class structures to connect metabolite/compound databases to empirical compounds;
and to connect experimental data to empirical compounds.

'''
import json
from operator import itemgetter
from mass2chem.epdsConstructor import epdsConstructor
from .search import *
from .ions import compute_adducts_formulae, generate_ion_signature

def read_table_to_peaks(infile, delimiter='\t'):
    '''Made as template. Modify accordign to your own file format.
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


#----------------------------------------------------------------------------------------
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

    def search_emp_cpd_single(self, emp_cpd, mode='pos', mz_tolerance_ppm=5):
        '''
        emp_cpd format has to follow (loosely) specificiations in metDataModel, requiring anchor ion, e.g.
        {interim_id': 12,
        'neutral_formula_mass': None,
        'neutral_formula': None,
        'MS1_pseudo_Spectra': [{'feature_id': 'F94',
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
            db_peaks = compute_adducts_formulae(R['neutral_formula_mass'], R['neutral_formula'], mode)
            results.append((_M['parent_epd_id'],
                            _M['ion_relation'],
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
        self.peak_to_empCpd = {}

    def build_from_list_peaks(self, list_peaks):
        '''
        list of peaks, e.g. [ {
                                'id_number': 555,
                                'mz': 133.0970, 
                                'apex': 654, 
                                'height': 14388.0, 
                                'left_base': 648, 
                                'right_base': 655, 
                                }, ...]
        '''
        self.list_peaks = list_peaks
        ECCON = epdsConstructor(list_peaks, mode=self.mode)

        list_empCpds = ECCON.peaks_to_epds()        # exclude_singletons=False)
        # updating epdsConstructor to handle singletons - done

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
        for interim_id, epd in self.dict_empCpds.items():
            peaks = epd['MS1_pseudo_Spectra']
            for P in peaks:
                self.peak_to_empCpd[P['feature_id']] = interim_id
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
            compound['neutral_formula_mass'], compound['neutral_formula'], self.mode)
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

    def annotate_against_KCD(self, KCD):
        '''
        Get all empCpd matches between this experimental dataset and a known compound database.
        KCD: knownCompoundDatabase instance.
        >>> KCD.search_emp_cpd_single( EED.dict_empCpds[15] )
        [('C6H14N_100.112624', 'M[1+]', 2), ('C6H13N_99.104799', 'M+H[1+]', 2)]

        '''
        resultDict = {}
        for epd in self.dict_empCpds.values():
            resultDict[epd['interim_id']] = KCD.search_emp_cpd_single(epd, self.mode)
        return resultDict

    def choose_top_epd(list_KCD_matches):
        '''recommend best from matches, e.g. 
        [('C6H14N_100.112624', 'M[1+]', 2), ('C6H13N_99.104799', 'M+H[1+]', 2)]
        by 1) highest score, and 2) peck order of ions (M+H > M+Na > M+ > M+H2O+H)
        '''
        pass

    def export_empCpds(self, outfile="EED_empCpds.json"):
        with open(outfile, 'w', encoding='utf-8') as O:
            json.dump(list(self.dict_empCpds.values()), O, ensure_ascii=False, indent=2)

    def export_annotations(self, search_result_dict, KCD, export_file_name_prefix, include_succinct=False):
        '''
        Export expt empCpds, match relationships, and a tsv table for input list_peaks, with all matched db empCpds.

        To-do:
        Will add option to include a version succinct annotation of anchor ions and recommended matches only.
        
        '''
        self.export_empCpds(export_file_name_prefix+"EED_empCpds.json")
        with open(export_file_name_prefix+"mapping.json", 'w', encoding='utf-8') as O:
            json.dump(search_result_dict, O, ensure_ascii=False, indent=2)

        s = "[peak]id_number\tmz\tapex\t[EmpCpd]interim_id\tmatched_DB_shorts\tmatched_DB_records\n"
        for peak in self.list_peaks:
            ii, mz, apex = peak['id_number'], peak['mz'], peak['apex']
            matched_DB_shorts, matched_DB_records = '', ''
            interim_id = self.peak_to_empCpd.get(ii, None)  # id on expt data; ii could be 0
            list_matches = []
            if interim_id != None:
                list_matches = search_result_dict[interim_id]
                # [('C7H8N2O_136.063663', 'M+H[1+]', 2), ('C7H9N2O_137.071488', 'M[1+]', 2), ('C7H6N2_118.053098', 'M+H2O+H[1+]', 1)]
            if not list_matches:
                list_matches = KCD.search_mz_single(mz)
                # [{'mz': 130.017306555, 'parent_epd_id': 'C4H3FN2O2_130.017856', 'ion_relation': 'M[1+]'}]
                list_matches = [(x['parent_epd_id'], x['ion_relation'], None) for x in list_matches]

            if list_matches:
                    matched_DB_shorts = ", ".join([ "(" + KCD.short_report_emp_cpd(xx[0]) + ")"  for xx in list_matches])
                    matched_DB_records = ", ".join([str(xx) for xx  in list_matches])

            s += '\t'.join([str(x) for x in [
                ii, mz, apex, interim_id, matched_DB_shorts, matched_DB_records]]) + "\n"

        outfile = export_file_name_prefix + '.tsv'
        with open(outfile, 'w') as O:
            O.write(s)

        print("\nAnnotation of %d Empirical compounds was written to %s." %(len(self.dict_empCpds), outfile))

