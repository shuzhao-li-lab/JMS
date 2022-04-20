
import pandas as pd
import yaml
import warnings
import json
from metDataModel.core import EmpiricalCompound

def read_metadata_file(input_file_path):
    '''
    Read a metadata YAML file.

    Input
    =====
    the input file path

    An example of a qualified YAML file
    ======================
    ion_relationship: '[M+H]+'
    mode: 'HILICpos'
    description: 'HILICpos authentic standard library from Shujian and Minghao'
    curator: 'Shujian Zheng|Minghao Gong'
    method: '5 min method'
    notes:
    initiation:
        curator: MG
        date: ''
        note: '5 min method, the library was created and combined by Minghao Gong and Shujian Zheng'

    version: '0.1'
    label: 'combn_1strnd'

    Output
    ======
    A dictionary of the metadata
    '''
    with open(input_file_path, 'r') as f:
        metadata = yaml.safe_load(f)
    return(metadata)
    
def libTab2EpdJson(df_path = '', delimiter = ', ',metadata = '', RT_long_format = True):
    '''
    Input
    =====
    df_path: the input tabular file path
    delimiter: the delimiter used in separate values in individual cells (note: this delimiter is not for separation between cells) 
    metadata: see also read_metadata_file
    RT_long_format: If there are incidents multiple retention time for a single compound stored in individual row

    Note
    ====
    the provided table needs to be comma separated data
    the metadata needs to be in yaml format
    The columns need to follow the rule as: 'mz','rt','identifier','name','Formula','adduct'
    If multiple adducts/mz of a particular empirical compound present, it needs to be at the same row (retention time is always the same with different adducts)
    The number of delimited items in identifier, name & Formula should be the same.

    Example of input tabular data
    =============================
    mz	rt	identifier	name	Formula	adduct
	130.0498696	231.24	HMDB0000805, HMDB0000267	D-Pyroglutamic acid, L-Pyroglutamic acid	C5H7NO3, C5H7NO3	[M+H]+
    130.0862551	127.38	HMDB0000716	L(-)-Pipecolinic acid	C6H11NO2	[M+H]+
    130.0862551	230.16	HMDB0000716	L(-)-Pipecolinic acid	C6H11NO2	[M+H]+
    132.0655197	154.14	HMDB0060460, HMDB0000725	cis-4-Hydroxy-D-proline, L-Hydroxyproline	C5H9NO3, C5H9NO3	[M+H]+
    132.0655197	219.3	HMDB0060460, HMDB0000725	cis-4-Hydroxy-D-proline, L-Hydroxyproline	C5H9NO3, C5H9NO3	[M+H]+

    Output
    ======
    a dictionary contains a list of empirical compound and a metadata

    '''
    df = pd.read_csv(df_path, dtype=str)

    if RT_long_format == True:
        reshaped_df = df.groupby([x for x in df.columns if 'rt'not in x]).agg({'rt': lambda x: ', '.join(str(i) for i in x)})
        reshaped_df = reshaped_df.reset_index()
    else:
        reshaped_df = df
    myEpds = []
    final_dict = reshaped_df.to_dict('index')
    
    for row in final_dict.values():
        Epd = EmpiricalCompound()
        
        # dealing with identity
        if delimiter in row['identifier']:
            ids = row['identifier'].split(delimiter)
            names = row['name'].split(delimiter)
            formulas = row['Formula'].split(delimiter)
        else:
            ids = [row['identifier']]
            names = [row['name']]
            formulas = [row['Formula']]
        for _id, name, formula in zip(ids,names,formulas):
            Epd.identity.append({'compounds': _id, 'names': name, 'formulas':formula, 'score': 0, 'probability': 0})
                   
        # dealing with MS1_pseudo_Spectra
        if delimiter in row['rt']:
            rts = row['rt'].split(delimiter)
        else:
            rts = [row['rt']]
        if delimiter in row['mz']:
            mzs = row['mz'].split(delimiter)
            adducts = row['adduct'].split(delimiter)
        else:
            mzs = [row['mz']]
            adducts = [row['adduct']]
        for rt in rts:
            for mz,adduct in zip(mzs,adducts):
                Epd.MS1_pseudo_Spectra.append({'feature_id': '', 
                                               'mz': float(mz.rstrip()), 
                                           'rtime': float(rt.rstrip()), 'charged_formula': '', 'ion_relation': adduct},)
        myEpds.append(Epd)
    
    Epd_dict = {}
    Epd_dict['list_of_Empirical_Compounds'] = [x.__dict__ for x in myEpds]
    Epd_dict['metadata'] = metadata

    return(Epd_dict)

def write_json(EpdJson_dict,output_file_path):
    '''
    write the empirical compound dictionary into output file under the correct path.
    '''
    with open(output_file_path, "w") as f:
        json.dump(EpdJson_dict,f,indent = 2)

def combn_lib(lib1,lib2,metadata = ''):
    '''
    Combine two library Json empirical compound dictionary 

    Input
    =====
    Two Empirical compound dictionary (the dictionary includes a list of empirical compound and a previous metadata)
    A metadata YAML file (See also read_metadata_file)

    Output
    ======
    Output a merged Empirical compound json file, attached with new metadata based on YAML file
    If there are overlaps of compound identifiers between two json files, create a warning message for further checking

    '''
    lib1_cpd_id = [cpd['compounds'] for epd in lib1['list_of_Empirical_Compounds'] for cpd in epd['identity']]
    lib2_cpd_id = [cpd['compounds'] for epd in lib2['list_of_Empirical_Compounds'] for cpd in epd['identity']]
    intersect = set(lib1_cpd_id).intersection(set(lib2_cpd_id))
    
    lib3 = {'list_of_Empirical_Compounds':  [], 'metadata' : metadata}
    if len(intersect) !=0:
        intersect_string = ";".join(list(intersect))
        warnings.warn(f"There are duplicated identifiers. Check for duplications of {intersect_string}")
    lib3['list_of_Empirical_Compounds'] = lib1['list_of_Empirical_Compounds'] + lib2['list_of_Empirical_Compounds']

    return(lib3)

   