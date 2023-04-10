'''
Conversion of genome scale metabolic models to JSON formats,and to indexed Python internal objects.
Matching btw metabolic models and experimental empirical compounds.
Used in mummichog 3.
'''

from khipu.epdsConstructor import epdsConstructor
from .dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

default_parameters = {'mode': 'pos',
                        'mz_tolerance_ppm': 5,
                        'rt_tolerance': 2,
                        'isotope_search_patterns': [(1.003355, '13C/12C', (0, 0.8)),
                                (2.00671, '13C/12C*2', (0, 0.8))],
                        'adduct_patterns': [(21.982, 'Na/H'),
                                (41.026549, 'ACN'),
                                (35.9767, 'HCl'),
                                (37.955882, 'K/H')],
                        'extended_adducts': [(1.0078, 'H'),
                                (-1.0078, '-H'),
                                (10.991, 'Na/H, double charged'),
                                (0.5017, '13C/12C, double charged'),
                                (117.02655, '-NH3'),
                                (17.02655, 'NH3'),
                                (-18.0106, '-H2O'),
                                (18.0106, 'H2O'),
                                (18.033823, 'NH4'),
                                (27.01089904, 'HCN'),
                                (37.94694, 'Ca/H2'),
                                (32.026215, 'MeOH'),
                                (43.96389, 'Na2/H2'),
                                (67.987424, 'NaCOOH'),
                                (83.961361, 'KCOOH'),
                                (97.96737927, 'H2SO4'),
                                (97.97689507, 'H3PO4')],
 }

def convert_json_model(jmodel):
    '''
    Convert a JSON style, concise genome scale metabolic models to mummichog style with internal indexing.
    All compound identifiers and other identifiers are expected to be contained within the model.
    Returns model as a dictionary with as in mummichog style.

    >>> jmodel.keys()
        dict_keys(['id', 'list_of_reactions', 'list_of_compounds', 'list_of_pathways', 'meta_data'])
    >>> jmodel['list_of_compounds'][772]
    {'id': 'ebastineoh', 'name': 'Hydroxylated Ebastine', 'identifiers': {'inchi': 'InChI=1S/C32H39NO3/c1-32(2,24-34)28-17-15-25(16-18-28)30(35)14-9-21-33-22-19-29(20-23-33)36-31(26-10-5-3-6-11-26)27-12-7-4-8-13-27/h3-8,10-13,15-18,29,31,34H,9,14,19-24H2,1-2H3'}, 'neutral_formula': '', 'charge': 0, 'charged_formula': 'C32H39NO3', 'neutral_mono_mass': 0.0, 'SMILES': '', 'inchi': ''}
    >>> jmodel['list_of_pathways'][2]
    {'id': 'group3', 'name': 'Aminoacyl-tRNA biosynthesis', 
    'list_of_reactions': ['HMR_5130', 'HMR_5131', 'HMR_5132', 'HMR_5133', 'HMR_5134', 'HMR_5135', 
            'HMR_5136', 'HMR_5137', 'HMR_5138', 'HMR_5139', 'HMR_5140', 'HMR_5141', 'HMR_5142', 'HMR_5143', 
            'HMR_5144', 'HMR_5145', 'HMR_5146', 'HMR_5147', 'HMR_5148', 'HMR_5149', 'HMR_5150']}
    >>> 
    >>> jmodel['list_of_reactions'][22]
    {'id': '24_25DHVITD3tm', 'reactants': ['2425dhvitd3'], 'products': ['2425dhvitd3']}
    >>> mfn['metabolic_pathways'][51]
    {'cpds': ['C00024', 'C00027', 'strdnccoa', 'C00100', 'C01944', 'q10h2', 'dmhptcoa', 'etfrd', 'C03035', 'etfox', 'dmnoncoa', 'arachdcoa', 'dcsptn1coa', 'ptdcacoa', 'C00412', 'od2coa', 'hdcoa', 'C02593', 'dlnlcgcoa', 'vacccoa', 'C00010', 'C00016', 'adrncoa', 'C00154', 'lneldccoa', 'lnlccoa', 'c226coa', 'lnlncacoa', 'C00399', 'odecoa', 'hpdcacoa', 'C01352', 'eicostetcoa', 'tmndnccoa', 'arachcoa'], 
    'rxns': ['FAOXC2242046m', 'FAOXC200180m', 'FAOXC204', 'FAOXC160', 'FAOXC16080m', 'FAOXC180', 'FAOXC80', 'FAOXC140', 'FAOXC1811601m', 'FAOXC1811603m', 'FAOXC150m', 'FAOXC18480m', 'FAOXC1811602m', 'FAOXC16180m', 'FAOXC170m', 'FAOXC18280m', 'FAOXC182806m', 'FAOXC183803m', 'FAOXC2252053m', 'FAOXC2031836m', 'FAOXC11', 'FAOXC204184m', 'FAOXC226205m', 'FAOXC226', 'ETFQO', 'ETF'], 
    'ecs': ['1.3.99.3', '1.5.5.1'], 'id': 'mfn1v10path136', 'name': 'Fatty acid oxidation'}
    >>> mfn['metabolic_rxns'][2]
    {'reactants': ['C06611'], 'id': 'R05233', 'source': 'kegg_dre', 'ecs': ['1.1.1.1'], 'products': ['C06613'], 'cpds': ['C06611', 'C06613'], 'pathway': '3-Chloroacrylic acid degradation'}
    >>>    
    mfn['cpd2pathways'] = {
         'C06125': ['Glycosphingolipid metabolism', 'Glycosphingolipid biosynthesis - ganglioseries'], ...
         }
    >>> mfn['cpd_edges'][:5]
[   ['C00106', 'C00526'], ['C05977', 'CE6326'], ['C00043', 'G00174'], ['C00199', 'C00231'], ['C00014', 'C00655']]
    '''
    new = {}
    new['id'] = jmodel['id']
    new['version'] = jmodel['meta_data']['version']
    cpdDict,  dict_cpds_def = {}, {}
    for cpd in jmodel['list_of_compounds']:
        if 'neutral_mono_mass' in cpd:
            cpd['neutral_formula_mass'] = cpd['neutral_mono_mass']
        if 'neutral_formula_mass' in cpd and cpd['neutral_formula_mass']:
            cpdDict[cpd['id']] = cpd    # convert_compound_azimuth_mummichog(cpd, hmdb_dict, kegg_dict)
            dict_cpds_def[cpd['id']] = cpd['name']
    new['Compounds'] = cpdDict
    new['dict_cpds_def'] = dict_cpds_def
    new['metabolic_rxns'] = jmodel['list_of_reactions']

    # bunch of indexing
    dict_rxns, edge2rxn, edge2enzyme = {}, {}, {}
    cpd_edges = []
    for rxn in jmodel['list_of_reactions']:
        dict_rxns[rxn['id']] = rxn
        _edges = []
        for cpd1 in rxn['reactants']:
            for cpd2 in rxn['products']:
                k = [cpd1, cpd2]
                cpd_edges.append(k)
                _edges.append(k)
                edge2rxn[','.join(sorted(k))] = rxn['id']
                enzymes = rxn.get('ecs', []) or rxn.get('enzymes', [])
                for e in enzymes:
                    edge2enzyme[','.join(sorted(k))] = e

    new['cpd_edges'] = cpd_edges
    new['edge2rxn'] = edge2rxn
    new['edge2enzyme'] = edge2enzyme
    # now do pathways
    metabolic_pathways = []
    for P in jmodel['list_of_pathways']:
        _p = {}
        _p['id'] = P['id']
        _p['name'] = P['name']
        _p['rxns'] = P['list_of_reactions']
        cpds, ecs, genes = [], [], []
        for ii in P['list_of_reactions']:
            rxn = dict_rxns[ii]
            cpds += rxn['reactants'] + rxn['products']
            enzymes = rxn.get('ecs', []) or rxn.get('enzymes', [])
            ecs += enzymes
            if 'genes' in rxn:
                genes += rxn['genes']
        _p['cpds'] = list(set(cpds))
        _p['ecs'] = list(set(ecs))
        _p['genes'] = list(set(genes))
        metabolic_pathways.append(_p)

    cpd2pathways = {}
    for _p in metabolic_pathways:
        for cpd in _p['cpds']:
            if cpd in cpd2pathways:
                cpd2pathways[cpd].append(_p['name'])
            else:
                cpd2pathways[cpd] = [_p['name']]

    new['metabolic_pathways'] = metabolic_pathways
    new['cpd2pathways'] = cpd2pathways

    return new



class metabolicNetwork:
    '''
    Simple wrapper class of a genome scale metabolic model in JSON.

    ?? keep or not ??
    '''
    def __init__(self, MetabolicModel):
        '''
        MetabolicModel : metabolic model in JSON style dictionary.
        '''
        #print_and_loginfo( "Loading metabolic network %s..." %MetabolicModel.version ) # version from metabolic model
        
        self.MetabolicModel = MetabolicModel
        self.network = self.build_network(MetabolicModel['cpd_edges'])
        
        self.version = MetabolicModel['version']
        self.Compounds = MetabolicModel['Compounds']
        self.metabolic_pathways = MetabolicModel['metabolic_pathways']
        self.dict_cpds_def = MetabolicModel['dict_cpds_def']
        self.cpd2pathways = MetabolicModel['cpd2pathways']
        self.edge2enzyme = MetabolicModel['edge2enzyme']
        self.total_cpd_list = self.network.nodes()
        
        
    def build_network(self, edges):
        '''
        import networkx as nx
        '''
        return nx.from_edgelist( edges )
        

    def get_pathways(self):
        pass

    def validate(self):
        pass


class DataMeetModel:
    '''
    Match the compounds in a genome scale metabolic model to a metabolomics experiment,
    via the format of empirical compounds. The empirical compounds are constructed via khipu.
    Their JSON format can accommodate added fields on the fly, and return results 
    as part of empCpd['identity'] (https://github.com/shuzhao-li/metDataModel).
    Because neutral mass is inferred from khipu, the match problem is simplified by 
    focusing on neutral mass (formula).
    '''
    def __init__(self, MetabolicModel, userFeatureList=None, userListEmpCpds=None, parameters=default_parameters):
        '''
        Besides MetabolicModel, this takes either userFeatureList or userListEmpCpds.
        parameters : dictionary to pass ion mode, m/z and rt tolerance and isotope/adduct patterns.
        MetabolicModel : metabolic model in JSON style dictionary.
        userFeatureList : list of JSON style features
        userListEmpCpds : list of empirical compounds, which can be constructed and processed elsewhere
        '''
        self.model = MetabolicModel
        self.userFeatureList = userFeatureList
        self.userListEmpCpds = userListEmpCpds
        
        self.mode = parameters['mode']
        self.isotope_search_patterns = parameters['isotope_search_patterns']
        self.adduct_patterns = parameters['adduct_patterns']
        self.extended_adducts = parameters['extended_adducts']
        self.mz_tolerance_ppm = parameters['mz_tolerance_ppm']
        self.rt_tolerance = parameters['rt_tolerance']

    def match_all(self):
        '''
        Match model compounds with empirical compounds, via the JMS KCD-EED architecture.
        The search here cannot distinguish isomers, thus centering on neutral mass/formula.

        Returns
        -------
        dict_empCpds : {id: empCpd, ...} with matched compounds in empCpd['list_matches'],
            which have KCD empCpd identifiers but compound records need to pull out KCD later.
            This dict includes singletons and empCpds without matches.

        Examples
        --------
        A intermediary EED.dict_empCpds before updating identity: 
        {'interim_id': 'kp203_202.1317', 
            'neutral_formula_mass': 202.13169603323, 'neutral_formula': None, 
            'Database_referred': [], 'identity': [], 
            'MS1_pseudo_Spectra': [{'id_number': 'F3900', 'mz': 204.1424, 'rtime': 23.04, 'rtime_left_base': '21.88', 
                    'rtime_right_base': '25.37', 'parent_masstrack_id': '2519', 'peak_area': '111965925', 'cSelectivity': '0.86', 
                    'goodness_fitting': '0.98', 'snr': '6675', 'detection_counts': '15', ' 'apex': 23.04, 
                    'representative_intensity': '111965925', 'id': 'F3900', 
                    'isotope': '13C/12C', 'modification': 'M+H+', 'ion_relation': '13C/12C,M+H+', 'parent_epd_id': 'kp203_202.1317'}, 
            {'id_number': 'F3684', 'mz': 203.1389, 'rtime': 23.04, 'rtime_left_base': '21.88', 'rtime_right_base': '25.37', 
                    'parent_masstrack_id': '2495', 'peak_area': '402792917', 'cSelectivity': '0.86', 'goodness_fitting': '0.98', 
                    'snr': '120', 'detection_counts': '15',  'apex': 23.04, 'representative_intensity': '402792917', 
                    'id': 'F3684', 'isotope': 'M0', 'modification': 'M+H+', 'ion_relation': 'M0,M+H+', 'parent_epd_id': 'kp203_202.1317'}], 
            'MS2_Spectra': [], 
            'list_matches': [('C9H18N2O3_202.131742', 'neutral', 1)]}

        An indexed KCD empCpd, KCD.mass_indexed_compounds['C9H18N2O3_202.131742']:
        {'interim_id': 'C9H18N2O3_202.131742',
            'neutral_formula': 'C9H18N2O3',
            'neutral_formula_mass': 202.13174244717,
            'compounds': [{'id': 'MAM03375',
            'name': 'L-Alanyl-L-Leucine',
            'identifiers': [['humanGEM', 'MAM03375'],
                ['bigg.metabolite', 'CE5866'],
                ['pubchem.compound', '6992388'],
                ['vmhmetabolite', 'CE5866'],
                ['metanetx.chemical', 'MNXM15786'],
                ['inchi',
                'InChI=1S/C9H18N2O3/c1-5(2)4-7(9(13)14)11-8(12)6(3)10/h5-7H,4,10H2,1-3H3,(H,11,12)(H,13,14)/t6-,7+/m1/s1']],
            'neutral_formula': 'C9H18N2O3',
            'charge': 0,
            'charged_formula': 'C9H18N2O3',
            'neutral_mono_mass': 202.13174244717,
            'SMILES': '',
            'inchi': '',
            'neutral_formula_mass': 202.13174244717}]}
        '''
        KCD = knownCompoundDatabase()
        KCD.mass_index_list_compounds(self.model['Compounds'].values())
        KCD.build_emp_cpds_index()
        EED = ExperimentalEcpdDatabase(mode=self.mode, 
                                       mz_tolerance_ppm=self.mz_tolerance_ppm, 
                                       rt_tolerance=self.rt_tolerance)
        if self.userFeatureList:
            EED.build_from_list_peaks(self.userFeatureList)
        elif self.userListEmpCpds:
            EED.build_from_list_empCpds(self.userListEmpCpds)
        EED.extend_empCpd_annotation(KCD)
        EED.annotate_singleton_mummichog(KCD)

        return self.update_identity(EED.dict_empCpds, KCD)


    def update_identity(self, dict_empCpds, KCD):
        '''
        Updates dict_empCpds by adding field `identity` with compound IDs.
        E.g. {'identity': ['MAM03375'], ... }
        This overwrites field `identity` if it exists.
        '''
        def _get_kcd_ids_(interim_id, KCD):
            e = KCD.mass_indexed_compounds.get(interim_id, None)
            if e:
                return [x['id'] for x in e['compounds']]
            else:
                return []

        for k,v in dict_empCpds.items():
            if 'list_matches' in v:
                v['identity'] = []
                for M in v['list_matches']:
                    v['identity'] += _get_kcd_ids_(M[0], KCD)

        return dict_empCpds


    def _construct_EmpiricalCompounds_(self):
        '''
        For testing only. Use EED class for full application, as in self.match_all().
        Returns
        -------
        dict_empCpds : {id: empCpd, ...}
        '''
        ECCON = epdsConstructor(self.userFeatureList, mode=self.mode)
        return ECCON.peaks_to_epdDict(
                        self.isotope_search_patterns,
                        self.adduct_patterns,
                        self.extended_adducts,
                        self.mz_tolerance_ppm,
                        self.rt_tolerance
        )
