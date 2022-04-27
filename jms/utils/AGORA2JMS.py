#!/usr/bin/env python
# coding: utf-8

# # Porting genome scale metabolic models for metabolomics
# 
# https://github.com/VirtualMetabolicHuman/AGORA/tree/master/CurrentVersion/AGORA_1_03
# 
# **Use cobra to parse SBML models whereas applicable**
# 
# **Base our code on metDataModel**
# 
# Each model needs a list of Reactions, list of Pathways, and a list of Compounds.
# It's important to include with Compounds with all linked identifiers to other DBs (HMDB, PubChem, etc), and with formulae (usually charged form in these models) when available.
# We can alwasy update the data later. E.g. the neural formulae can be inferred from charged formula or retrieved from public metabolite database (e.g., HMDB) if linked.
# Save in Python pickle and in JSON.
# 
# **No compartmentalization**
# - After decompartmentalization,
#   - transport reactions can be removed - they are identified by reactants and products being the same.
#   - redundant reactions can be merge - same reactions in diff compartments become one.
# Georgi Kolishovski, Minghao Gong, 2022-04-21

# !pip install cobra --user --ignore-installed ruamel.yaml
# !pip install --upgrade metDataModel # https://github.com/shuzhao-li/metDataModel/ 
# !pip install --upgrade numpy pandas

import cobra # https://cobrapy.readthedocs.io/en/latest/io.html#SBML
from metDataModel.core import Compound, Reaction, Pathway, MetabolicModel
import pandas as pd
import sys
import os


sys.path.append("/Users/gongm/Documents/projects/mass2chem/")
sys.path.append("/Users/gongm/Documents/projects/JMS/JMS/JMS")
from mass2chem.formula import *
from jms.formula import *
from jms.utils.gems import *
from jms.utils.git_download import git_download_from_direcotry

from datetime import datetime
today =  str(datetime.today()).split(" ")[0]

def port_metabolite(M):
    '''convert cobra Metabolite to metDataModel Compound'''
    Cpd = Compound()
    Cpd.src_id = remove_compartment_by_split(M.id,'[') # remove the [c] from eg h2o[c]
    Cpd.id = remove_compartment_by_split(M.id,'[') # remove the [c] from eg h2o[c]
    Cpd.name = M.name
    Cpd.charge = M.charge
    Cpd.charged_formula = M.formula
    Cpd.neutral_formula = adjust_charge_in_formula(M.formula,M.charge)
    Cpd.neutral_mono_mass = neutral_formula2mass(Cpd.neutral_formula)
    Cpd.db_ids = list(M.notes.items())
    mydict = M.notes   # other databaseIDs  are in the notes tag
    Cpd.SMILES= mydict.get("SMILES",None) # not know if this is useful or not
    Cpd.inchi= mydict.get("InChIString",None)
    return Cpd

# port reactions
def port_reaction(R):
    '''port reaction'''
    new = Reaction()
    new.id = R.id
    new.reactants = [remove_compartment_by_split(m.id,'[') for m in R.reactants] 
    new.products = [remove_compartment_by_split(m.id,'[') for m in R.products] 
    return new

# Port pathway
def port_pathway(P):
    new = Pathway()
    new.id = P.id
    new.source = ['AGORA',]
    new.name = P.name
    new.list_of_reactions = [x.id for x in P.members]
    return new

if __name__ == '__main__':
    # The path you intended to save
    input_dir = '/Users/gongm/Documents/projects/JMS/JMS/JMS/test/input/test_automatic_AGORA' # 
    output_fdr = '/Users/gongm/Documents/projects/JMS/JMS/JMS/test/input/test_automatic_AGORA' #
    version_manual = 'AGORA_1_03_With_Mucins_sbml' #

    download_from_git = False
    if download_from_git == True:
        path_in_github = 'CurrentVersion/AGORA_1_03/AGORA_1_03_With_Mucins_sbml' #

    # download the AGORA xml files from github
    if download_from_git == True:
        user = 'VirtualMetabolicHuman'
        repo_name = 'AGORA'
        git_download_from_direcotry(user,repo_name,path_in_github,input_dir)
    
    # load the xml files
    file_names = [x for x in os.listdir(input_dir) if x.endswith('.xml')]
    for file_name in file_names:

        try:
            version = path_in_github.split("/")[-1]
        except:
            version = version_manual

        note = f'AGORA cloned from https://github.com/VirtualMetabolicHuman, retrieved from {today}\ .'

        metadata = {
                    'species': file_name.split(".")[0],
                    'version': version,
                    'sources': [f'AGORA, retrieved {today}'], #
                    'status': '',
                    'last_update': today,  #
                    'note': note
                }

        file_path = os.path.join(input_dir, file_name)
        model = cobra.io.read_sbml_model(file_path)
        myCpds = []
        for i in range(len(model.metabolites)):
            myCpds.append(port_metabolite(model.metabolites[i]))

        print(f'Before decompartmentalization, there are {len(myCpds)} compounds')

        # remove duplicated compounds
        myCpds = remove_duplicate_cpd(myCpds)

        print(f'After decompartmentalization, there are {len(myCpds)} compounds left')


        myCpds = fetch_AGORA_GEM_identifiers(myCpds,json_path = '../data/staged/vmh.json',overwrite = True)

        ## Reactions to port
        myRxns = []
        for R in model.reactions:
            myRxns.append( port_reaction(R) )

        print(f'Before removing transport reactions, there are {len(myRxns)} reactions')

        # remove duplicated reactions after decompartmentalization
        myRxns = remove_duplicate_rxn(myRxns)

        print(f'After removing transport reactions, there are {len(myRxns)} reactions')


        ## Pathways to port
        myPathways = []
        for P in model.groups:
            myPathways.append(port_pathway(P))

        # retain the valid reactions in list of pathway
        myPathways = retain_valid_Rxns_in_Pathways(myPathways,myRxns)

        print(f'There are {len(myPathways)} groups or pathways in the model')

        # Collected data; now output

        ## metabolicModel to export
        MM = MetabolicModel()
        MM.meta_data = metadata
        species = metadata['species']
        MM.id = f'{species}' 
        MM.list_of_pathways = [P.serialize() for P in myPathways]
        MM.list_of_reactions = [R.serialize() for R in  myRxns]
        MM.list_of_compounds = [C.serialize() for C in myCpds]

        # Write pickle file
        export_pickle(os.path.join(output_fdr,f'{MM.id}.pickle'), MM)
        print(f'Export pickle file')

        # Write json file
        export_json(os.path.join(output_fdr,f'{MM.id}.json'), MM)
        print(f'Export json file')

        # Write dataframe 

        export_table(os.path.join(output_fdr,f'{MM.id}_list_of_compounds.csv'),MM, 'list_of_compounds')
        print(f'Export a table of the list of compounds')


