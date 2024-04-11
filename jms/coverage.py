'''
Test coverage of a dataset or a model on pathways/ontologies.

Could exclude currency metabolites in future.

'''

from .modelConvert import convert_json_model, DataMeetModel
from .search import build_centurion_tree, find_all_matches_centurion_indexed_list
from .data.humangem_pathways import humangem
# from .data.smpdb_pathways import smpdb
from .empiricalCpds import *



def report_khipu_statistics(json_epd_file, snr=5, shape=0.9, natural_ratio_limit=0.5):
    '''
    Get numbers of khipus and singletons, and isopairs and good khipus.
    Good khipu = with isopair and M0 being a good feature.

    num_isopair_mtracks is based on num_good_khipus. 
    Thus the ratio num_good_khipus/num_isopair_mtracks indicates ave good isomers observed.
    
    This is for data without isotope labeling. 
    Planing to use function report_khipu_statistics_from_labeled_data for labeled data.
    
    Returns stats_dict, good_khipus.

    '''
    list_epds = list(json.load(open(json_epd_file)).values())
    num_features, num_good_features = 0, 0
    for epd in list_epds:
        for f in epd["MS1_pseudo_Spectra"]:
            num_features += 1
            f['snr'] = float(f['snr'])
            f['goodness_fitting'] = float(f['goodness_fitting'])
            f['peak_area'] = float(f['peak_area'])
            f['is_good_peak'] = check_good_peak(f, snr, shape)
            if f['is_good_peak']:
                num_good_features += 1
    
    khipus_isopairs, num_isopair_mtracks, good_khipus = get_isopairs_good_khipus(list_epds, natural_ratio_limit)
    num_empcpds = len(list_epds)
    num_singletons = count_singletons(list_epds)
    num_good_khipus = len(good_khipus)
    if num_isopair_mtracks:
        ratio_good_khipus_mtracks = num_good_khipus/num_isopair_mtracks
    else:
        ratio_good_khipus_mtracks = None
    
    stats_dict = {
        'num_features': num_features,
        'num_good_features': num_good_features,
        'num_empcpds': num_empcpds,
        'num_khipus': num_empcpds - num_singletons,
        'num_singletons': num_singletons,
        'num_khipus_isopairs': len(khipus_isopairs),
        'num_good_khipus': num_good_khipus,
        'num_isopair_mtracks': num_isopair_mtracks,
        'ratio_good_khipus_mtracks': ratio_good_khipus_mtracks
    }
    
    return stats_dict, good_khipus


def report_khipu_statistics_from_labeled_data(isotope_labeling_ratio=1.2):
    '''
    
    This analysis is often sample dependent, as control samples are usually included in the expt.
    
    see .empiricalCpds.filter_khipus_by_samples
    
    '''
    pass



def get_feature_stats_per_table(infile, snr=5, shape=0.9, sep='\t'):
    '''
    Get num_features, num_good_features from asari feature table.
    Alternatively, one can get stats from the empCpd JSON file.
    '''
    def _is_good_feature(line, snr, shape, sep):
        a = line.split(sep)
        _snr, _shape = float(a[7]), float(a[6])
        if _snr > snr and _shape > shape:
            return True
        else:
            return False
        
    lines = open(infile).readlines()
    num_features = len(lines) - 3
    num_good_features = len([line for line in lines[3:] if _is_good_feature(line, snr, shape, sep)])

    return num_features, num_good_features
    

def check_good_peak(peak, snr=5, shape=0.9, area=0):
    '''
    Decide if peak is good quality based on SNR and peak shape.
    Not using peak area by default value.
    '''
    if peak['snr'] > snr and peak['goodness_fitting'] > shape and peak['peak_area'] > area:
        return True
    else:
        return False


def report_pathway_coverage_from_gem(model, list_neutral_epds, outfile):
    '''
    This writes a report outfile and returns metabolic_pathways
    One can choose to use all neutral epds or only good epds.
    e.g. 
    stats_dict, good_khipus = report_khipu_statistics(json_epd_file, snr=5, shape=0.9, natural_ratio_limit=0.5)
    neutrals = get_neutrals( good_khipus )
    
    model is the JSON GEM model, e.g.
    model = json.load(open('metabolicModel_az_HumanGEM_20220302_noCompartmentalization.json'))
    model.keys()    
    dict_keys(['id', 'list_of_reactions', 'list_of_compounds', 'list_of_pathways', 'meta_data'])

    report_pathway_coverage modifies the pathways - thus pathways are reloaded here everytime.
    '''
    mcgmodel = convert_json_model(model)
    metabolic_pathways = report_pathway_coverage_from_neutrals( mcgmodel, list_neutral_epds)
    export_pathway_coverage_table(metabolic_pathways, outfile)
    return metabolic_pathways


def report_pathway_coverage_from_neutrals(mcg_model, list_neutral_epds, ppm=5):
    '''Returns a list of pathways with matched cpds to list_epds.

    list_epds : list of empirical compounds, but mz is neutral mass.
    mcg_model : metabolic model in mummichog style.

    '''
    mztree = build_centurion_tree(list_neutral_epds)
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



def get_pathwayCoverage_cpds_to_epds(list_epds, pathway_collection=humangem, ppm=5):
    '''
    Match list_epds to list_of_metabolites in full metabolic model (i.e. pathway_collection).
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
    This test each pathway in mcg_model and collect matched epds.

    list_epds : list of empirical compounds.
    mcg_model : metabolic model in mummichog style.

    returns metabolic_pathways.
    '''
    neutrals = get_neutrals( 
                            filter_epds(list_epds, neutral_formula_mass, multiple_ions) 
                            )
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
    '''
    There are functions handling empCpd filtering in report_pathway_coverage.
    They are preferred to be done explicitly outside this.
    See report_pathway_coverage_from_gem
    '''
    list_epds = load_epds_from_json(epd_json_file)
    mcg_model = convert_json_model(json_model)
    metabolic_pathways = report_pathway_coverage(
        mcg_model, list_epds, neutral_formula_mass, multiple_ions, ppm)
    export_pathway_coverage_table(metabolic_pathways, outfile)




#
# ------------ test ---------------
#

if __name__ == "__main__":
    
    model = json.load(open('/Users/lish/li.proj/ipsc_atlas/metabolicModel_az_HumanGEM_20220302_noCompartmentalization.json'))
      
    t = '/Users/lish/li.workspace/MT_AML_13C/mzML_HILICneg_20240318_THP1_MOLM13_asari_project_32693656/Annotated_empiricalCompounds.json'
    stats_dict, good_khipus = report_khipu_statistics(t, snr=5, shape=0.9, natural_ratio_limit=0.5)
    neutrals = get_neutrals( good_khipus )
    print("\n\n~~~~~~~~~ num good neutrals ~~~~~~~~~~\n\n")
    print(len(neutrals))
    print(stats_dict)
    
    report_pathway_coverage_from_gem(model, neutrals, outfile='test_jms_coverage.tsv')
    