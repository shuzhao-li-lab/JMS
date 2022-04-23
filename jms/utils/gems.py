import sys
import pickle
import json
import pandas as pd
from jms.formula import *

sys.path.append("/Users/gongm/Documents/projects/mass2chem/")
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

def neutral_formula2mass(neutral_formula):
    '''
    Convert neutral formula to mass but removing characters (e.g., X, R) typical in GEM but not ready for metabolomics application
    '''
    formula_dict = parse_chemformula_dict(neutral_formula)
    if ("R" not in formula_dict) & ("X" not in formula_dict) & (len(formula_dict) != 0) :
        mono_mass = calculate_mass(formula_dict,6)
    else:
        mono_mass = None
    return(mono_mass)

def remove_compartment(Cpd_id,len_of_suffix):
    '''
    remove compartment based the length of suffix of the compound ID.
    '''
    Cpd_id_mod = Cpd_id.rstrip()[:len(Cpd_id)-len_of_suffix]
    return(Cpd_id_mod)

def remove_duplicate_cpd(list_of_Cpds):
    '''
    remove duplicated compounds after removing the compartment information.
    '''
    id_list = []
    new_list_of_Cpds = []
    for Cpd in list_of_Cpds:
        if Cpd.id not in id_list:
            id_list.append(Cpd.id)
            new_list_of_Cpds.append(Cpd)
        else:
            continue
    return(new_list_of_Cpds)

def remove_duplicate_rxn(list_of_Rxns):
    '''
    Removing uniport or exchange reactions.
    '''
    new_list_of_Rxns = []
    for Rxn in list_of_Rxns:
        if Rxn.reactants == Rxn.products:
            continue
        else:
            new_list_of_Rxns.append(Rxn)
    return(new_list_of_Rxns)

def retain_valid_Rxns_in_Pathways(list_of_pathways, list_of_reactions):
    '''
    retain only valid Reactions in the pathway list after transport reactions excluded.
    '''
    valid_Rxns_list = [x.id for x in list_of_reactions]
    new_list_of_pathways = []
    for path in list_of_pathways:
        new_list_of_reactions = []
        for rxn in path.list_of_reactions:
            if rxn in valid_Rxns_list:
                new_list_of_reactions.append(rxn)
        path.list_of_reactions = new_list_of_reactions
        new_list_of_pathways.append(path)
    return(new_list_of_pathways)

def export_pickle(export_file_path,MetabolicModel):
    with open(export_file_path, 'wb') as f:
        pickle.dump(MetabolicModel.serialize(), f, pickle.HIGHEST_PROTOCOL)

def export_json(export_file_path,MetabolicModel):
    s = json.JSONEncoder().encode( MetabolicModel.serialize() )
    with open(export_file_path, 'w') as f:
        f.write(s)

def export_table(export_file_path = '',MetabolicModel = '',list_of_entries = 'list_of_compounds'):
    pd.DataFrame(MetabolicModel.serialize()[list_of_entries]).to_csv(export_file_path, index = False)

