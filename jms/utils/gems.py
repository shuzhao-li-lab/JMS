import sys
import pickle
import json
import pandas as pd
import re

from jms.utils.git_download import * 
from jms.formula import *

from mass2chem.formula import *

def list_all_identifiers(cobra_model,notes):
    '''
    List all the identifier sources and the example of identifiers in the current model
    It is used to make strategy for parsing. 

    Input
    =====
    cobra_model: a cobra model
    notes: dependent on the identifier source, it may differ. In AGORA and EBI models, 'notes' are by default while in humanGEM, 'annotation' are by default

    Output
    ======
    A set of unique identifier source and examples of the identifiers

    '''
    identifiers = []
    example = {}
    for item in cobra_model.metabolites:
        try:
            identifiers.extend(item.__dict__[notes].keys())
            for k,v in item.__dict__[notes].items(): 
                if k not in example:
                    example.update({k:v})
        except:
            None
    uni_id = set(identifiers)
    return(uni_id,example)


def dict2listOfTuple(db_dict, delimiter):
    '''
    Convert dictionary to a list of tuple for convenient conversion in db_ids
    '''
    list_of_tuple = []
    for k,v in db_dict.items():
        if delimiter in v:
            list_of_ids = [x.rstrip() for x in v.split(delimiter)]
            for _id in list_of_ids:
                list_of_tuple.append((k,_id))
        else:
            list_of_tuple.append((k,v))
    return list_of_tuple

def listOfTuple2dict(tuples, delimiter = None):
    '''
    Convert a list of tuple to dictionary for convenient conversion in db_ids
    '''
    res_dict = {}
    unique_db_name = list(set([tup[0] for tup in tuples]))
    for db_name in unique_db_name:
        db_ids = []
        for tup in tuples:
            if tup[0] == db_name:
                db_ids.append(tup[1])
        if len(db_ids) == 1:
            db_ids = db_ids[0]
            res_dict[db_name] = db_ids
        elif len(db_ids) > 1:
            if delimiter == None:
                res_dict[db_name] = db_ids
            else:
                res_dict[db_name] = ', '.join(db_ids)
    return res_dict

def neutral_formula2mass(neutral_formula:str)->float:
    """Convert neutral formula to mass but removing characters (e.g., X, R) typical in GEM but not ready for metabolomics application

    Parameters
    ----------
    neutral_formula : str
        Neutral Formula

    Returns
    -------
    float
        Mass of the neutral formula
    Examples
    --------
    >>> from jms.utils.gems import *
    >>> neutral_formula2mass('C17H32O2')
    268.24023
     
    """
    formula_dict = parse_chemformula_dict(neutral_formula)
    if ("R" not in formula_dict) & ("X" not in formula_dict) & ("Z" not in formula_dict) & (len(formula_dict) != 0) :
        mono_mass = calculate_mass(formula_dict,6)
    else:
        mono_mass = float(0)
    return mono_mass

def remove_compartment_by_substr(Cpd_id:str,len_of_suffix:int)->str:
    """    
    remove compartment based the length of suffix of the compound ID.
    This is working especially for humanGEM as they just use the last character for compartmentalization designation
    However, this is not working for EBI model as they do have "_bm" alone with "_c"

    Parameters
    ----------
    Cpd_id : str
        Compound id from GEMs db
    len_of_suffix : int
        Length of suffix need to cut

    Returns
    -------
    str
        Compound id after remove

    Examples
    --------
    >>> remove_compartment_by_substr('MAM00001c', 1)
    'MAM00001'  

    """
    Cpd_id_mod = Cpd_id.rstrip()[:len(Cpd_id)-len_of_suffix]
    return Cpd_id_mod

def remove_compartment_by_split(Cpd_id,delimiter):
    '''
    remove compartment based the delimiters of the compound ID.
    This works for most other models.
    '''
    Cpd_id_mod = Cpd_id.rsplit(delimiter,1)[0]
    return(Cpd_id_mod)

def remove_duplicate_cpd(list_of_Cpds:list)->list:
    """remove duplicated compounds after removing the compartment information.

    Parameters
    ----------
    list_of_Cpds : list
        Raw compound list

    Returns
    -------
    list
        Compound list after deduplication
     
    """
    id_list = []
    new_list_of_Cpds = []
    for Cpd in list_of_Cpds:
        if Cpd.id not in id_list:
            id_list.append(Cpd.id)
            new_list_of_Cpds.append(Cpd)
        else:
            continue
    return new_list_of_Cpds

def remove_duplicate_rxn(list_of_Rxns:list)->list:
    """Removing uniport or exchange reactions.

    Parameters
    ----------
    list_of_Rxns : list
        Raw reaction list

    Returns
    -------
    list
        Reaction list after deduplication
     
    """
    new_list_of_Rxns = []
    for Rxn in list_of_Rxns:
        if Rxn.reactants == Rxn.products:
            continue
        else:
            new_list_of_Rxns.append(Rxn)
    return new_list_of_Rxns

def retain_valid_Rxns_in_Pathways(list_of_pathways:list, list_of_reactions:list)->list:
    """retain only valid Reactions in the pathway list after transport reactions excluded.

    Parameters
    ----------
    list_of_pathways : list
        list of pathways after port_pathway()
    list_of_reactions : list
        list of reactions after port_reaction()

    Returns
    -------
    list
        list of pathways after transport reactions excluded
     
    """
    valid_Rxns_list = [x.id for x in list_of_reactions]
    new_list_of_pathways = []
    for path in list_of_pathways:
        new_list_of_reactions = []
        for rxn in path.list_of_reactions:
            if rxn in valid_Rxns_list:
                new_list_of_reactions.append(rxn)
        path.list_of_reactions = new_list_of_reactions
        new_list_of_pathways.append(path)
    return new_list_of_pathways

def export_pickle(export_file_path:str,MetabolicModel):
    """Export metDataModel.core.MetabolicModel to pickle

    Parameters
    ----------
    export_file_path : str
        pickle path
    MetabolicModel : metDataModel.core.MetabolicModel
        metabolic model to export
     
    """
    with open(export_file_path, 'wb') as f:
        pickle.dump(MetabolicModel.serialize(), f, pickle.HIGHEST_PROTOCOL)

def export_json(export_file_path:str,MetabolicModel):
    """Export metDataModel.core.MetabolicModel to json

    Parameters
    ----------
    export_file_path : str
        json path
    MetabolicModel : metDataModel.core.MetabolicModel
        metabolic model to export
     
    """
    s = json.JSONEncoder().encode( MetabolicModel.serialize() )
    with open(export_file_path, 'w') as f:
        f.write(s)

def export_table(export_file_path:str,MetabolicModel,list_of_entries:str):
    """Export compounds, reactions and pathways in metabolic model into .csv

    Parameters
    ----------
    export_file_path : str
        csv file path
    MetabolicModel : metDataModel.core.MetabolicModel
        the model store metabolites, reactions and pathways
    list_of_entries : str
        suffix of .csv file
     
    """
    pd.DataFrame(MetabolicModel.serialize()[list_of_entries]).to_csv(export_file_path, index = False)

def fetch_MetabAtlas_GEM_identifiers(compound_list,
                                     modelName,
                                     local_path,
                                     metab_file_name = 'metabolites.tsv',
                                     overwrite = True):
    url_met = f'https://github.com/SysBioChalmers/{modelName}/blob/main/model/metabolites.tsv'
    git_download_from_file(url_met,local_path,metab_file_name)
    metab_df = pd.read_csv(os.path.join(local_path,metab_file_name),sep = '\t')
    metab_df.index = metab_df.mets
    metab_dict = metab_df.to_dict('index')

    metab_id_list = list(set([x[0:len(x)-1] for x in metab_dict.keys()]))
    decomp_metab_dict = {}
    for Id_decomp in metab_id_list:
        for Id_comp,v in metab_dict.items():
            if Id_decomp in Id_comp:
                decomp_metab_dict.update({Id_decomp:v})
                break
    for Cpd in compound_list:
        for k,v in decomp_metab_dict.items():
            new_db_ids = []
            if Cpd.id == k:
                if overwrite == True:
                    for kk,vv in v.items():
                        if re.match('met(.*)ID',kk) and isinstance(vv,str):
                            sub_key = re.match('met(.*)ID',kk)[1]
                            if ';' in vv:
                                sub_list = vv.split(';')
                                for sub_id in sub_list:
                                    new_db_ids.append((sub_key,sub_id))
                            else:
                                new_db_ids.append((sub_key,vv))
                    Cpd.db_ids = new_db_ids

def fetch_AGORA_GEM_identifiers(compound_list,
                                json_path,
                                overwrite = True):
    with open(json_path,'r') as f:
        list_vmh_cpd = json.load(f)
    vmh_dict = {}
    for vmh_cpd in list_vmh_cpd:
        vmh_dict.update({vmh_cpd['id']:vmh_cpd})
    new_cpd_list = []
    for myCpd in compound_list:
        for k,v in vmh_dict.items():
            if myCpd.id == k:
                if overwrite == True:
                    myCpd.db_ids = v['identifiers'] # the vmh json is using `identifiers` rather than `db_ids`
                break
        new_cpd_list.append(myCpd)
    return new_cpd_list


def fetch_CarveMe_GEM_charge_formula(compound_list,Bigg_json_path,overwrite = True):
    with open(Bigg_json_path,'r') as f:
        list_Bigg_cpd = json.load(f)
    Bigg_dict = {}
    for Bigg_cpd in list_Bigg_cpd:
        Bigg_dict.update({Bigg_cpd['id']:Bigg_cpd})
    new_cpd_list = []
    for myCpd in compound_list:
        for k,v in Bigg_dict.items():
            if myCpd.id == k:
                if overwrite == True:
                    myCpd.charge = v['charge']
                    myCpd.charged_formula = v['charged_formula'] # the Bigg json is using `identifiers` rather than `db_ids`
                break
        new_cpd_list.append(myCpd)
    return new_cpd_list