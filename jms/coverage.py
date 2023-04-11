'''
Test coverage of a dataset or a model on pathways/ontologies.

Could exclude currency metabolites in future.
Refactoring...

'''

from .modelConvert import convert_json_model, DataMeetModel
from .search import build_centurion_tree, find_all_matches_centurion_indexed_list
from .data.humangem_pathways import humangem
# from .data.smpdb_pathways import smpdb
from .empiricalCpds import *

def get_pathwayCoverage_cpds_to_epds(list_epds, pathway_collection=humangem, ppm=5):
    '''
    Match list_epds to list_of_metabolites with pathway memberships.
    list_epds : input list of empirical compounds. Users can apply filtering before this.
    pathway_collection : humangem = {
                                'version': 'https://github.com/SysBioChalmers/Human-GEM, retrieved 2022-02-09',
                                'dict_pathways': dict_pathways,
                                'list_of_metabolites': list_of_metabolites,
                            }
    Returns
    -------
    List of matches as [(cpd['id'], cpd['name'], matched_list of epds), ...]

    Examples
    --------
    list_of_metabolites[55]:
        {'id': 'MAM00067', 'mw': 1101.43877624496, 'name': '(2E)-tricosenoyl-CoA', 'pathways': ['Fatty acid elongation (odd-chain)']}
    list(humangem['dict_pathways'].items())[11]:
        ('group12',
        {'name': 'Beta oxidation of branched-chain fatty acids (mitochondrial)',
        'cpds': ['MAM00845',
        'MAM01802',
        ...,
        'MAM02040']})
    '''
    neutrals = get_neutrals(list_epds)
    mztree = build_centurion_tree(neutrals)
    matched = []
    for cpd in pathway_collection['list_of_metabolites']:
        MM = find_all_matches_centurion_indexed_list(cpd['mw'], mztree, ppm)
        if MM:
            matched.append((cpd['id'], cpd['name'], MM))

    return matched

def get_pathwayCoverage_mummichog(list_epds, metabolic_model, ppm=5):
    '''Match list_epds to pathways in mummichog style, using modelConvert.DataMeetModel.
    '''
    DMM = DataMeetModel(MetabolicModel=metabolic_model, userFeatureList=None, userListEmpCpds=list_epds)
    DMM.mz_tolerance_ppm = ppm
    return DMM.match_all()


#
# ------------------------------------------------------------------------------------------
# Convenience functions below for quick analysis; not suitabe for production software 

def report_pathway_coverage(mcg_model, list_epds, neutral_formula_mass=True, multiple_ions=True, ppm=5):
    '''Returns a list of pathways with matched cpds to list_epds.

    list_epds : list of empirical compounds.
    mcg_model : metabolic model in mummichog style.

    '''
    neutrals = get_neutrals( filter_epds(list_epds, neutral_formula_mass, multiple_ions) )
    mztree = build_centurion_tree(neutrals)
    for p in mcg_model['metabolic_pathways']:
        cpds = [mcg_model['Compounds'].get(x, None) for x in p['cpds']]
        cpds = [x for x in cpds if x]
        for c in cpds:
            c['mw'] = c.get('neutral_mono_mass', None) or c.get('neutral_formula_mass', None)
        
        MM = get_match(cpds, mztree, ppm)
        p['number_cpds'] = len(p['cpds'])   # likely greater than cpds
        p['number_matched_epds'] = len(MM)
        p['matched_epds'] = [(x['id'], x['name'], x['mw'], y) for x,y in MM]
        p['cpd_names'] = [mcg_model['dict_cpds_def'].get(x, x) for x in p['cpds']]

    return mcg_model['metabolic_pathways']


def export_pathway_coverage_table(metabolic_pathways, outfile='pathway_matches.tsv'):
    '''Export stats table from metabolic_pathways, after running coverage analysis.
    '''
    s = 'pid\tname\tN_matched\tN_in_pathway\tcpds_in_pathway\tmatched_epds\n'
    for p in metabolic_pathways:
        s += '\t'.join([p['id'], p['name'], str(p['number_matched_epds']), str(p['number_cpds']),
                        str(p['matched_epds']), 
                        str(p['cpd_names']), 
                        ]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)


def wrapper_file2coverage(json_model, epd_json_file, 
                          neutral_formula_mass=True, multiple_ions=True, ppm=5,
                          outfile='pathway_matches.tsv'):
    list_epds = load_epds_from_json(epd_json_file)
    mcg_model = convert_json_model(json_model)
    metabolic_pathways = report_pathway_coverage(
        mcg_model, list_epds, neutral_formula_mass, multiple_ions, ppm)
    export_pathway_coverage_table(metabolic_pathways, outfile)
