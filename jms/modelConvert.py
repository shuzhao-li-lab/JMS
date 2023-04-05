'''
Conversion of genome scale metabolic models to JSON formats,and to indexed Python internal objects.
Matching btw metabolic models and experimental empirical compounds.
Used in mummichog 3.
'''

from khipu.epdsConstructor import epdsConstructor
from .dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase


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
    def __init__(self, parameters, MetabolicModel, userFeatureList):
        '''
        parameters : dictionary to pass ion mode, m/z and rt tolerance and isotope/adduct patterns.
        MetabolicModel : metabolic model in JSON style dictionary.
        userFeatureList : list of JSON style features
        '''
        self.model = MetabolicModel
        self.userFeatureList = userFeatureList
        
        self.mode = parameters['mode']
        self.isotope_search_patterns = parameters['isotope_search_patterns']
        self.adduct_patterns = parameters['adduct_patterns']
        self.extended_adducts = parameters['extended_adducts']
        self.mz_tolerance_ppm = parameters['mz_tolerance_ppm']
        self.rt_tolerance = parameters['rt_tolerance']

    def _get_ListOfEmpiricalCompounds_(self):
        '''
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

    def match_all(self):
        '''
        Match model compounds with empirical compounds, via the JMS KCD-EED architecture.
        The search here cannot distinguish isomers, thus centering on neutral mass/formula.

        Returns
        -------
        dict_empCpds : {id: empCpd, ...} with matched compounds in empCpd['list_matches']
        '''
        KCD = knownCompoundDatabase()
        KCD.mass_index_list_compounds(self.model['Compounds'])
        KCD.build_emp_cpds_index()
        EED = ExperimentalEcpdDatabase(mode=self.mode, 
                                       mz_tolerance_ppm=self.mz_tolerance_ppm, 
                                       rt_tolerance=self.rt_tolerance)
        EED.dict_empCpds = self._get_ListOfEmpiricalCompounds_()
        EED.extend_empCpd_annotation(KCD)
        EED.annotate_singletons(KCD)       
        return EED.dict_empCpds

