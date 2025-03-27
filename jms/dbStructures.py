'''
Class structures to connect metabolite/compound databases to empirical compounds;
and to connect experimental data to empirical compounds.
'''

import json
from operator import itemgetter
import numpy as np
import tqdm
from khipu.epdsConstructor import epdsConstructor

from khipu.utils import adduct_search_patterns, \
                            adduct_search_patterns_neg, \
                                isotope_search_patterns, \
                                    extended_adducts

from .search import *
from .formula import *
from .ions import compute_adducts_formulae, generate_ion_signature
from .data.list_formula_mass import list_formula_mass


def annotate_peaks_against_kcds(list_peaks, 
                                list_compounds, 
                                export_file_name_prefix='jms_annotated_',
                                mode='pos',  
                                mz_tolerance_ppm=5, 
                                rt_tolerance = 2,
                                check_isotope_ratio = False):
    '''
    Wrapper function as example, to generate three annotation files for input list_peaks.
    list_compounds is known compound database, e.g.
    list_compounds = json.load(open('jms/data/compounds/list_compounds_HMDB4.json'))
    list_peaks is experimental data, list of features in json.
    '''
    KCD = knownCompoundDatabase()
    KCD.mass_index_list_compounds(list_compounds)
    KCD.build_emp_cpds_index()
    # optional to export_mass_indexed_compounds
    # KCD.export_mass_indexed_compounds(export_file_name_prefix+"KCD_mass_indexed_compounds.json")

    EED = ExperimentalEcpdDatabase(mode=mode, 
                                   mz_tolerance_ppm=mz_tolerance_ppm,
                                   rt_tolerance=rt_tolerance)
    EED.build_from_list_peaks(list_peaks)
    # Second, singletons that get a formula match in KCD
    EED.annotate_singletons(KCD)       
    # Third, the remaining features unmatched to anything (orphans). Exported for potential downstream work.
    # EED.dict_empCpds = self.append_orphans_to_epmCpds(EED.dict_empCpds)
    
    EED.export_annotations(KCD, export_file_name_prefix)


#----------------------------------------------------------------------------------------
class knownCompoundDatabase:
    '''
    Indexed data store for known compounds, e.g., from a database or a metabolic model. 
    The search should be on empCpd not Cpd.
    An empCpd is an empirical compound, a set of features drived from the same mass (see README), 
    which often include isomers. The regular mass search cannot distinguish isomers.

    One can search by mass or mass tree (patterns of isotopes/adducts, in the form of empCpd.
    We use build_centurion_tree to index empCpds;
    centurion_mass_tree is an indexed dictionary to group ions by 100th decimal.
    There are situations where only neutral mass is searched; others requiring ionized forms.
    This class generates three trees to accommodate them: 'neutral', 'pos' and 'neg'.
    '''
    def __init__(self):
        '''
        Two main data containers, mass_indexed_compounds and emp_cpds_trees.
        The latter is indexed for searches, separately for positive and negative ion modes.
        '''
        self.mass_indexed_compounds = {}
        self.emp_cpds_trees = { 'pos': {}, 
                                'neg': {},
                                'neutral': {},
                                }
        
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
        from collections import defaultdict
        _db = defaultdict(list)
        for cpd in list_compounds:
            _db[cpd['neutral_formula'] + "_" + str(round(float(cpd['neutral_formula_mass']),6))].append(cpd)
        for k, v in _db.items():
            self.mass_indexed_compounds[k] = {
                "interim_id": k,
                "neutral_formula": v[0]['neutral_formula'],
                "neutral_formula_mass": float(v[0]['neutral_formula_mass']),
                "compounds": v
            }

    def build_emp_cpds_index(self, primary_only=True, include_C13=True):
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
        
        peak_lists = {"pos": [], "neg": [], "neutral": []}
        for k, v in tqdm.tqdm(self.mass_indexed_compounds.items()):
            peak = {a: b for a, b in v.items()}
            peak['mz'] = v['neutral_formula_mass']
            peak['parent_epd_id'] = k
            peak['ion_relation'] = 'neutral'
            peak_lists['neutral'].append(peak)
            for mode in ["pos", "neg"]:
                for ion in __ion_generator__(v['neutral_formula_mass'], v['neutral_formula'], mode=mode, primary_only=primary_only):
                    ion_peak = dict(peak)
                    ion_peak['mz'] = ion[0]
                    ion_peak['ion_relation'] = ion[1]
                    if "order" in ion:
                        ion_peak["order"]: ion[3]
                    peak_lists[mode].append(ion_peak)
        self.emp_cpds_trees = {k: build_centurion_tree(v) for k, v in peak_lists.items()}

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
        emp_cpd : empirical compound, following loosely specificiations in metDataModel, e.g.
        {interim_id': 12,
        'neutral_formula_mass': None,
        'neutral_formula': None,
        'MS1_pseudo_Spectra': [{'id_number': 'F94',
                                'mz': 98.97531935309317,
                                'rtime': 700.0,
                                'charged_formula': '',
                                'ion_relation': 'anchor',
                                'parent_epd_id': 12}, ], ...}

        If input emp_cpd is generated from khipu, neutral_formula_mass is usually alrady assigned, 
        and this is simple search on neutral_formula_mass.
        Otherwise, anchor ion is searched in emp_cpds_trees (from mass_indexed_compounds).
        The ion_relations are expected to be ['M[1+]', 'M+H[1+]', 'M+Na[1+]', 'M+H2O+H[1+]'] for pos,
        and ['M[-]', 'M-H[-]', 'M-H2O-H[-]', 'M+Cl[-]'] for neg.

        This searches only neutral mass or anchor ion, because if anchor ion is not found in the database, 
        other adducts are not possible.

        Return
        ======
        list of matched db empCpd as [(empCpd_id, ion_relation, score), ...], sorted by score.
        While most search will return a single match. Occasionally, undeterministic situations appear,
        e.g. C6H13N and C6H14N could match to M+ and M+H+, respectively.
        The resolution of these will be left to downstream methods (bayesian etc).
        The score is number of matched ions, which may or may not be useful.
        '''
        results = []
        if emp_cpd['neutral_formula_mass']:
            matches = self.search_mz_single(emp_cpd['neutral_formula_mass'], 'neutral', mz_tolerance_ppm)
            for _M in matches:
                results.append((_M['parent_epd_id'],
                                _M['ion_relation'], 1))
        else:
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
    Build a boutique data store for user's experimental data, list of JSON features/peaks.
    Allow search of Compounds, which support targeted compound search in a dataset.
    '''
    def __init__(self, 
                 mode='pos', 
                 mz_tolerance_ppm=5, 
                 rt_tolerance=2,
                 ):
        '''
        mode: ionizaation mode, 'pos' or 'neg'.
        mz_tolerance_ppm: ppm tolerance in examining m/z patterns.
        rt_tolerance: tolerance threshold for deviation in retetion time, arbitrary unit depending on input data.
                Default intended as 2 seconds.
        The isotope/adduct patterns here are populated via self.get_isotope_adduct_patterns, as imported from khipu.
        But one can directly overwrite the patterns for customization (e.g. in asari).
        '''
        self.mode = mode
        self.mz_tolerance_ppm = mz_tolerance_ppm
        self.rt_tolerance = rt_tolerance
        self.get_isotope_adduct_patterns()
        self.list_peaks = []
        self.indexed_peaks = {}
        self.dict_peaks = {}
        self.dict_empCpds = {}
        self.indexed_empCpds = {}
        self.peak_to_empCpd = {}
        self.peak_to_empCpd_ion_relation = {}

    def get_isotope_adduct_patterns(self, 
                                    adduct_search_patterns=adduct_search_patterns,
                                    adduct_search_patterns_neg=adduct_search_patterns_neg,
                                    isotope_search_patterns=isotope_search_patterns[:2],
                                    extended_adducts=extended_adducts
                                    ):
        '''
        Populate isotope/adduct patterns for this instance.
        The default isotope/adduct patterns here are imported from khipu.
        The default isotope_search_patterns are limited to M0 and M1.
        '''
        self.adduct_patterns = adduct_search_patterns
        if self.mode == 'neg':
            self.adduct_patterns = adduct_search_patterns_neg
        self.isotope_search_patterns = isotope_search_patterns
        self.extended_adducts = extended_adducts

    def build_from_list_peaks(self, 
                              list_peaks, 
                              charges=[1, 2, 3],
                              has_parent_masstrack=True,
                              ):
        '''
        Wrapper of khipu epdsConstructor, to construct empirical compounds (empCpds).
        Updates self.dict_empCpds and does self.index_empCpds().

        list_peaks : [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'rtime': 654, 
                'height': 14388.0, 'id_number': 555}, ...]
        isotope/adduct patterns, mz_tolerance_ppm and rt_tolerance are from self,
        as they are parameters to this class.

        charges : multiple charged ions to consider.
        has_parent_masstrack : True if from asari. 
            Used to determine the need of setting parent_masstrack.

        Updates
        -------
        self.list_peaks
        self.dict_empCpds
        self.indexed_peaks, self.formula_tree, self.indexed_empCpds

        Future updates can include enforcing overlap in scans not just rtime tolerance.
        '''
        # check if intensity value is provided
        if 'representative_intensity' not in list_peaks[0]:
            if 'peak_area' in list_peaks[0]:
                for p in list_peaks:
                    p['representative_intensity'] = p['peak_area']
            elif 'height' in list_peaks[0]:
                for p in list_peaks:
                    p['representative_intensity'] = p['height']
            else:
                print("Warning: cannot get representative_intensity..")
                for p in list_peaks:
                    p['representative_intensity'] = 0

        self.list_peaks = list_peaks
        ECCON = epdsConstructor(list_peaks, mode=self.mode)
        self.dict_empCpds = ECCON.peaks_to_epdDict(
                        self.isotope_search_patterns,
                        self.adduct_patterns,
                        self.extended_adducts,
                        self.mz_tolerance_ppm,
                        self.rt_tolerance,
                        charges,
                        has_parent_masstrack,
        ) 
        self.index_empCpds()

    def build_from_list_empCpds(self, list_empCpds):
        for E in list_empCpds:
            self.dict_empCpds[E['interim_id']] = E
        self.index_empCpds()

    def index_empCpds(self):
        '''
        Build indices for self.list_peaks and self.empCpds.
        Khipu has peak:empCpd built N:1 but other methods should comply too.

        Updates
        -------
        self.indexed_peaks, self.formula_tree, self.indexed_empCpds
        '''
        for P in self.list_peaks:
            self.dict_peaks[P['id_number']] = P

        self.indexed_peaks  = build_centurion_tree(self.list_peaks)
        # use formula.get_formula_ions_tree
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

    def annotate_empCpds_against_KCD(self, KCD, mz_tolerance_ppm=5):
        '''
        Get all empCpd matches between this experimental dataset and a known compound database.
        KCD: knownCompoundDatabase instance. KCD.search_emp_cpd_single() searches only neutral mass or anchor ion. 
        Returns dictionary e.g.
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
            resultDict[epd['interim_id']] = KCD.search_emp_cpd_single(epd, self.mode, mz_tolerance_ppm )
        return resultDict

    def search_mz_for_formula(self, mz, mz_tolerance_ppm=5):
        '''
        return best matched formula or None, 
        format as {'mz': ion[0], 'neutral_formula': formula, 
                               'ion_relation': ion[1],
                               'neutral_formula_mass': mass}
        '''
        return search_mz_formula_tree(mz, self.formula_tree, limit_ppm=mz_tolerance_ppm)

    # Annotation functions
    def extend_empCpd_annotation(self, KCD):
        '''
        This searches KCD for the list of empCpds.
        With khipu, most empCpds have neutral mass, which is used for KCD search.
        If no neutral mass, an anchor ion is used for KCD search (KCD.search_emp_cpd_single).
        Returned list_matches have KCD empCpd identifiers but compound records need to pull out KCD later.
        '''
        epd_search_result_dict = self.annotate_empCpds_against_KCD(KCD, self.mz_tolerance_ppm)
        for interim_id, V in epd_search_result_dict.items():
            if V:
                self.dict_empCpds[interim_id]['list_matches'] = V

    def empCpds_formula_search(self, KCD):
        '''
        Update self.dict_empCpds for formulae first by KCD search then .data.formula_tree.
        KCD: knownCompoundDatabase instance.
        '''
        epd_search_result_dict = self.annotate_empCpds_against_KCD(KCD, self.mz_tolerance_ppm)
        for interim_id, V in epd_search_result_dict.items():
            if V:
                # update empCpd with formula matches
                self.dict_empCpds[interim_id]['list_matches'] = V
                #       expecting probem from khipu to this
                self.dict_empCpds[interim_id]['neutral_formula'] = \
                    KCD.mass_indexed_compounds[V[0][0]]['neutral_formula']
                self.dict_empCpds[interim_id]['neutral_formula_mass'] = \
                    KCD.mass_indexed_compounds[V[0][0]]['neutral_formula_mass']
            else:
                # first peak should be anchor
                anchor_mz = self.dict_empCpds[interim_id]['MS1_pseudo_Spectra'][0]['mz']
                F = self.search_mz_for_formula(anchor_mz, self.mz_tolerance_ppm)
                if F:
                    self.dict_empCpds[interim_id]['neutral_formula'] = F['neutral_formula']
                    self.dict_empCpds[interim_id]['neutral_formula_mass'] = F['neutral_formula_mass']

    def __extend_empCpds__(self):
        '''
        Extend empCpds by 
        additional isotopes/adducts, based on mass2chem.formula.compute_adducts_formulae

        Hold off; no need now from khipu results.
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
                    matches = find_all_matches_centurion_indexed_list(ion[0], peakTree, self.mz_tolerance_ppm)
                    for peak in matches:
                        if is_coeluted(anchor, peak, rt_tolerance=self.rt_tolerance):
                            peak['ion_relation'] = ion[1]
                            E['MS1_pseudo_Spectra'].append(peak)
                            self.peak_to_empCpd_ion_relation[peak['id_number']] = peak['ion_relation']
                            self.peak_to_empCpd[peak['id_number']] = E['interim_id']
                            # this may overwrite 1:N relationships

    def singleton_formula_search(self, KCD):
        '''
        Search singletons for formulae first by KCD search then .data.formula_tree.
        KCD: knownCompoundDatabase instance.

        Returns
        -------
        List of first matched KCD compound for each singleton in self.dict_peaks, e.g.
        [('F17904', {'interim_id': 15321, 'neutral_formula_mass': 916.606449, 'neutral_formula': 'C56H84O10',...}), ...]
        '''
        found = []
        singletons = [p for p in self.dict_peaks if p not in self.peak_to_empCpd.keys()]
        for p in singletons:
            _mz = self.dict_peaks[p]['mz']
            list_matches = KCD.search_mz_single(_mz, self.mode, self.mz_tolerance_ppm)
            # [{'mz': 130.017306555, 'parent_epd_id': 'C4H3FN2O2_130.017856', 'ion_relation': 'M[1+]'}]
            if list_matches:
                # take 1st match only here; will model better in future
                _epd = KCD.mass_indexed_compounds[list_matches[0]['parent_epd_id']]
                _epd['isotope'] = '13C/12C'
                _epd['ion_relation'] = _epd['modification'] = list_matches[0]['ion_relation']
                found.append((p, _epd))
            else:
                # formula search
                F = self.search_mz_for_formula(_mz, self.mz_tolerance_ppm)
                if F:
                    found.append((p, F))            # p is id_number

        return found
    

    def annotate_singleton_mummichog(self, KCD):
        '''
        This applies to KCD using a metabolic model as the sole database.

        Updates
        -------
        self.dict_empCpds : {id: empCpd, ...}
        '''
        singletons = [p for p in self.dict_peaks if p not in self.peak_to_empCpd.keys()]
        for p in singletons:
            peak = self.dict_peaks[p]
            interim_id = 'epd_' + peak['id_number']
            list_matches = KCD.search_mz_single(peak['mz'], self.mode, self.mz_tolerance_ppm)
            # [{'mz': 130.017306555, 'parent_epd_id': 'C4H3FN2O2_130.017856', 'ion_relation': 'M[1+]'}]
            if list_matches:
                self.dict_empCpds[interim_id] = {'interim_id': interim_id,
                    'MS1_pseudo_Spectra': [peak],
                    'list_matches': [(LL['parent_epd_id'], LL['ion_relation'], 1) for LL in list_matches],
                }


    def annotate_singletons(self, KCD):
        '''
        Search singletons for formulae first by KCD search then .data.formula_tree.
        Add new empCpds to self.dict_empCpds as new empCpds, and extend adduct search.
        Assign new identifiers as (10000 +number of dict_empCpds + numerical count).

        Updates
        -------
        self.dict_empCpds : {id: empCpd, ...}
        '''
        formula_to_peaks = {}
        found = self.singleton_formula_search(KCD)
        # multiple singletons can match to same formula as different ions, 
        # if ion calculation differs btw here and KCD
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

        new_id_start = len(self.dict_empCpds) + 10000
        for formula, PP in formula_to_peaks.items():
            neutral_formula_mass = PP[0][1]['neutral_formula_mass']
            P1 = self.dict_peaks[PP[0][0]]
            tmp = [P1, ]
            for jj in range(1, len(PP)):
                _P = self.dict_peaks[PP[jj][0]]
                if is_coeluted(tmp[-1], _P, rt_tolerance=self.rt_tolerance):
                    tmp.append(_P)
                else:
                    # not coeluted, new empCpd
                    new_id_start += 1
                    self.dict_empCpds[new_id_start] = {'interim_id': str(new_id_start),
                            'neutral_formula_mass': neutral_formula_mass, 'neutral_formula': formula,
                            'MS1_pseudo_Spectra': self.__extend_peakList__(
                                formula, neutral_formula_mass, tmp, peakTree, self.mz_tolerance_ppm),
                    }
                    tmp = [_P, ]

            new_id_start += 1
            self.dict_empCpds[new_id_start] = {'interim_id': str(new_id_start),
                    'neutral_formula_mass': neutral_formula_mass, 'neutral_formula': formula,
                    'MS1_pseudo_Spectra': self.__extend_peakList__(
                                formula, neutral_formula_mass, tmp, peakTree, self.mz_tolerance_ppm),
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
                if is_coeluted(anchor, peak, rt_tolerance=self.rt_tolerance):
                    peak['ion_relation'] = ion[1]
                    new.append(peak)
        return epd_peaks + new


    def annotate_all_against_KCD(self, KCD, mz_tolerance_ppm=5):
        '''
        Not used now.

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
                list_matches = KCD.search_mz_single(mz, self.mode, mz_tolerance_ppm)
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


    def export_annotations(self, KCD, export_file_name_prefix):
        '''
        Export expt empCpds and a tsv table for input list_peaks, with all matched db empCpds.
        '''
        self.export_empCpds(export_file_name_prefix+"EED_empCpds.json")
        s = "[peak]id_number\tmz\trtime\tapex(scan number)\t[EmpCpd]interim_id\t[EmpCpd]ion_relation\tneutral_formula\tneutral_formula_mass\
        \tname_1st_guess\tmatched_DB_shorts\tmatched_DB_records\n"
        
        for _, V in self.dict_empCpds.items():
            name_1st_guess, matched_DB_shorts, matched_DB_records = '', '', ''
            if 'list_matches' in V:
                list_matches = V['list_matches']
                if list_matches:
                    name_1st_guess = KCD.mass_indexed_compounds[list_matches[0][0]]['compounds'][0]['name']
                    matched_DB_shorts = ", ".join([ "(" + KCD.short_report_emp_cpd(xx[0]) + ")"  for xx in list_matches])
                    matched_DB_records = ", ".join([str(xx) for xx  in list_matches])

            for peak in V['MS1_pseudo_Spectra']:
                s += '\t'.join([str(x) for x in [
                    peak['id_number'], peak['mz'], peak['rtime'], peak['apex'], V['interim_id'], peak.get('ion_relation', ''),
                    V['neutral_formula'], V['neutral_formula_mass'],
                    name_1st_guess, matched_DB_shorts, matched_DB_records]]) + "\n"

        outfile = export_file_name_prefix + '.tsv'
        with open(outfile, encoding='utf-8', mode='w') as O:
            O.write(s)

        print("\nAnnotation of %d Empirical compounds was written to %s." %(len(self.dict_empCpds), outfile))



    #-----------------------------------------

    def update_annotation(self, KCD, peak_result_dict, adduct_patterns, mz_tolerance_ppm):
        '''
        Placeholder. Not used now but for future consieration -

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


            self.dict_empCpds.values(), remaining_peaks, adduct_patterns, mz_tolerance_ppm=mz_tolerance_ppm,
        ) 
        # update self.dict_empCpds
        for EC in list_empCpds:
            self.dict_empCpds[EC['interim_id']] = EC
