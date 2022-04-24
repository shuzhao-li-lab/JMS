#!/usr/bin/env python
# coding: utf-8

# # Porting genome scale metabolic models for metabolomics
# 
# **Human-GEM as default human model, for better compatibility**
# https://github.com/SysBioChalmers/Human-GEM
# 
# **Use cobra to parse SBML models whereas applicable**
# 
# Not all models comply with the formats in cobra. Models from USCD and Thiele labs should comply.
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
# 
# Shuzhao Li, 2021-10-21|
# Minghao Gong, 2022-04-19

# !pip install cobra --user --ignore-installed ruamel.yaml
# !pip install --upgrade metDataModel # https://github.com/shuzhao-li/metDataModel/ 
# !pip install --upgrade numpy pandas

import cobra # https://cobrapy.readthedocs.io/en/latest/io.html#SBML
from metDataModel.core import Compound, Reaction, Pathway, MetabolicModel
import requests
import sys
import os
from datetime import datetime

sys.path.append("/Users/gongm/Documents/projects/mass2chem/")
sys.path.append("/Users/gongm/Documents/projects/JMS/JMS/JMS")
from mass2chem.formula import *
from jms.formula import *
from jms.utils.gems import *

# The path you intended to save
output_fdr = '../../testdata/HumanGEM/'
name_xml = 'Human-GEM.xml'
github_xml_path = 'https://github.com/SysBioChalmers/Human-GEM/blob/main/model/Human-GEM.xml'

# metadata
today =  str(datetime.today()).split(" ")[0]

note = """Human-GEM compartmentalized, with genes and ECs."""

metadata = {
            'species': 'human',
            'version': '',
            'sources': [f'https://github.com/SysBioChalmers/Human-GEM, retrieved {today}'], #
            'status': '',
            'last_update': today,  #
            'note': note,
        }



# download the most updated Human-GEM.xml
HG_xml_path = os.path.join(output_fdr, name_xml)
with open(HG_xml_path, 'w') as f:
    r = requests.get(f'{github_xml_path}?raw=true')
    f.write(r.text)

# Read the model via cobra
xmlFile = HG_xml_path
model = cobra.io.read_sbml_model(xmlFile)

# Port metabolite

def port_metabolite(M):
    # convert cobra Metabolite to metDataModel Compound
    Cpd = Compound()
    Cpd.src_id = remove_compartment_by_substr(M.id,1)
    Cpd.id = remove_compartment_by_substr(M.id,1)              # temporarily the same with the source id
    Cpd.name = M.name
    Cpd.charge = M.charge
    Cpd.neutral_formula = adjust_charge_in_formula(M.formula,M.charge)
    Cpd.neutral_mono_mass = neutral_formula2mass(Cpd.neutral_formula)
    Cpd.charged_formula = M.formula
    Cpd.db_ids = [['humanGEM',Cpd.src_id]] # using src_id to also reference humanGEM ID in db_ids field
    for k,v in M.annotation.items():
        if k != 'sbo':
            if isinstance(v,list):
                Cpd.db_ids.append([[k,x] for x in v])
            else: 
                if ":" in v:
                    Cpd.db_ids.append([k,v.split(":")[1]])
                else:
                    Cpd.db_ids.append([k,v])
    
    inchi_list = [x[1].split('=')[1] for x in Cpd.db_ids if x[0] == 'inchi']
    if len(inchi_list) ==1:
        Cpd.inchi = inchi_list[0]
    elif len(inchi_list) >1:
        Cpd.inchi = inchi_list

    return Cpd

myCpds = []
for i in range(len(model.metabolites)):
    myCpds.append(port_metabolite(model.metabolites[i]))

print(f'Before decompartmentalization, there are {len(myCpds)} compounds')

# remove duplicated compounds
myCpds = remove_duplicate_cpd(myCpds)

print(f'After decompartmentalization, there are {len(myCpds)} compounds left')


# Port reactions
# port reactions, to include genes and enzymes
def port_reaction(R):
    new = Reaction()
    new.id = R.id
    new.reactants = [remove_compartment_by_substr(m.id,1) for m in R.reactants] # decompartmentalization
    new.products = [remove_compartment_by_substr(m.id,1) for m in R.products]   # decompartmentalization
    new.genes = [g.id for g in R.genes]
    ecs = R.annotation.get('ec-code', [])
    if isinstance(ecs, list):
        new.enzymes = ecs
    else:
        new.enzymes = [ecs]       # this version of human-GEM may have it as string
    return new

## Reactions to port
myRxns = []
for R in model.reactions:
    myRxns.append( port_reaction(R) )

print(f'Before removing transport reactions, there are {len(myRxns)} reactions')

# remove duplicated reactions after decompartmentalization
myRxns = remove_duplicate_rxn(myRxns)

print(f'After removing transport reactions, there are {len(myRxns)} reactions')

# ## Port pathway
# pathways, using group as pathway. Other models may use subsystem etc.

def port_pathway(P):
    new = Pathway()
    new.id = P.id
    new.source = ['Human-GEM v1.10.0',]
    new.name = P.name
    new.list_of_reactions = [x.id for x in P.members]
    return new

## Pathways to port
myPathways = []
for P in model.groups:
    myPathways.append(port_pathway(P))

print(f'There are {len(myPathways)} groups or pathways in the model')


# retain the valid reactions in list of pathway
myPathways = retain_valid_Rxns_in_Pathways(myPathways,myRxns)

# Collected data; now output

## metabolicModel to export
MM = MetabolicModel()
MM.id = f'az_HumanGEM_{today}' #
MM.meta_data = metadata
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
import pandas as pd
export_table(os.path.join(output_fdr,f'{MM.id}_list_of_compounds.csv'),MM, 'list_of_compounds')
export_table(os.path.join(output_fdr,f'{MM.id}_list_of_reactions.csv'),MM, 'list_of_reactions')
export_table(os.path.join(output_fdr,f'{MM.id}_list_of_pathways.csv'),MM, 'list_of_pathways')

print(f'Export a table of the list of compounds')

