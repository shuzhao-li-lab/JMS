'''
Conversion of genome scale metabolic models to JSON formats,and to indexed Python internal objects.

Started in mummichog 2.6, but moving to this repo and will continue with mummichog 3.

'''


import json
import numpy as np
import networkx as nx

from .ions import compute_isotopes_adducts

from .JSON_metabolicModels import metabolicModels
from .compound_dicts import hmdb_dict, kegg_dict

# 
# for porting GEMs
# see https://github.com/shuzhao-li/Azimuth/tree/master/docs
#



def json_convert_azmuth_mummichog(jmodel):
    '''
    The jmodel was parsed from one of the Genome scale metabolic models. 
    All compound identifiers and other identifiers are expected to be contained within the model.
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
        cpdDict[cpd['id']] = convert_compound_azimuth_mummichog(cpd, hmdb_dict, kegg_dict)
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
                # enzymes not in current test model and not tested ------------------ watch out for this, could be rxn['ecs']
                if 'enzymes' in rxn:
                    for e in rxn['enzymes']:
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
        cpds, ecs = [], []
        for ii in P['list_of_reactions']:
            rxn = dict_rxns[ii]
            cpds += rxn['reactants'] + rxn['products']
            # enzymes not in current test model and not tested ------------------ watch out for this, could be rxn['ecs']
            if 'enzymes' in rxn:
                ecs += rxn['enzymes']
        _p['cpds'] = list(set(cpds))
        _p['ecs'] = list(set(ecs))
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


def convert_compound_azimuth_mummichog(cpd, hmdb_dict, kegg_dict):
    '''
    In mummichog 2, compound takes the format of
    cpd format: {'formula': 'C29H51N2O8PRS', 'mw': 0, 'name': 'linoelaidic acid ACP (all trans)', 'adducts': {}}
    So mw = neutral_mono_mass; formula = neutral_formula.
    We built compound dictionaries on HMDB and KEGG separately, as id: (formula, mass).
    Return compound in mummichog 2 format.
    '''
    new = {}
    new['id'] = cpd['id']
    new['name'] = cpd['name']
    new['formula'] = cpd['neutral_formula']
    new['mw'] = cpd['neutral_mono_mass']
    if new['formula'] and new['mw'] > 1:        # both valid
        return new
    else:
        hmdb = cpd['identifiers'].get('hmdb', '')
        if hmdb:
            try:
                new['formula'], new['mw'] = hmdb_dict[hmdb]
            except KeyError:
                pass
        else:
            kegg = cpd['identifiers'].get('kegg.compound', '')
            if kegg:
                try:
                    new['formula'], new['mw'] = kegg_dict[kegg]
                except KeyError:
                    pass
        return new
