"""
Function used to port CUT models into metDataModel.

Yuanye Chi 2023-03-23
"""
import sys
import os
from datetime import datetime

import cobra  # https://cobrapy.readthedocs.io/en/latest/io.html#SBML

import requests

import logging

from mass2chem.formula import *

from jms.formula import *
from jms.port.port_utils import *
from jms.port.port_config import *

from metDataModel.core import Compound, Reaction, Pathway, MetabolicModel

# Logger Configuration
formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logging.basicConfig(level=logging.INFO, handlers=[handler])


def port_CUT(species: str):
    """port different species in CUT model

    Parameters
    ----------
    species : str
        species in CUT model. e.g. human, mouse, rat, yeast, worm...

    """
    model_source = Sources.CUT
    cur_info = basic_info_dict[Sources.CUT][species]

    logging.info(f""""
                ##################################
                #### start porting {species}  ####
                ##################################  
    """)

    # download the most updated XXX-GEM.xml
    download_path = os.path.join(cur_info['output_dir'], cur_info['name_xml'])

    logging.info(cur_info['github_xml_path'] + '?raw=true')
    with open(download_path, 'w') as f:
        r = requests.get(cur_info['github_xml_path'] + '?raw=true')
        f.write(r.text)

    # read xml by cobra
    model = cobra.io.read_sbml_model(download_path)

    ##################################
    #####     port metabolites   #####
    ##################################
    myCpds = [port_metabolite(metabolite, model_source)
              for metabolite in model.metabolites]
    logging.info(
        f'Before decompartmentalization, there are {len(myCpds)} compounds')
    # remove duplicated compounds
    myCpds = remove_duplicate_cpd(myCpds)
    logging.info(
        f'After decompartmentalization, there are {len(myCpds)} compounds left')

    ##################################
    #####     port reactions     #####
    ##################################
    myRxns = [port_reaction(reaction, model_source)
              for reaction in model.reactions]
    logging.info(
        f'Before removing transport reactions, there are {len(myRxns)} reactions')
    # remove duplicated reactions after decompartmentalization
    myRxns = remove_duplicate_rxn(myRxns)
    logging.info(
        f'After removing transport reactions, there are {len(myRxns)} reactions')

    ##################################
    #####     port pathways      #####
    ##################################
    # Pathways to port
    myPathways = [port_pathway(pathway, model_source, species)
                  for pathway in model.groups]
    logging.info(
        f'There are {len(myPathways)} groups or pathways in the model')
    # retain the valid reactions in list of pathway
    myPathways = retain_valid_Rxns_in_Pathways(myPathways, myRxns)

    ##################################
    #####     build model        #####
    ##################################
    # metabolicModel to export
    MM = MetabolicModel()
    MM.id = f'{species}_GEM_{today}'
    MM.meta_data = cur_info['metadata']
    MM.list_of_pathways = [P.serialize() for P in myPathways]
    MM.list_of_reactions = [R.serialize() for R in myRxns]
    MM.list_of_compounds = [C.serialize() for C in myCpds]

    ##################################
    #####      output model      #####
    ##################################
    # Write pickle file
    export_pickle(os.path.join(cur_info['output_dir'], f'{MM.id}.pickle'), MM)
    logging.info(f'Export pickle file')

    # Write json file
    export_json(os.path.join(cur_info['output_dir'], f'{MM.id}.json'), MM)
    logging.info(f'Export json file')

    # Write dataframe
    import pandas as pd
    export_table(os.path.join(
        cur_info['output_dir'], f'{MM.id}_list_of_compounds.csv'), MM, 'list_of_compounds')
    export_table(os.path.join(
        cur_info['output_dir'], f'{MM.id}_list_of_reactions.csv'), MM, 'list_of_reactions')
    export_table(os.path.join(
        cur_info['output_dir'], f'{MM.id}_list_of_pathways.csv'), MM, 'list_of_pathways')

    logging.info(f'Export a table of the list of compounds')

    logging.info(f""""
                ##################################
                #### finish porting {species} ####
                ##################################  
    """)
