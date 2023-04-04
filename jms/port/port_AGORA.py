"""
Function used to port AGORA models into metDataModel.

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

def git_download_from_direcotry(dir_url:str,local_output_dir:str)->list:
    """Download files from a github directory

    Parameters
    ----------
    dir_url : str
        the path of the github directory
    local_output_dir : str
        the path of the local directory you will to download

    Returns
    -------
    list
        list of paths of each xml file
     
    """
    list_of_attrs =json.loads(requests.get(dir_url).text)
    list_of_dicts = []
    for attrs in list_of_attrs:
        try:
            list_of_dicts.append(
                {'url':attrs['html_url'],
                'name': attrs['name']
                }
            ) # don't use url, which will not work as it go through git API and file size < 1 MB
        except:
            None

    logging.info('%d species found in AGORA', len(list_of_dicts))

    # download the urls
    for dict in list_of_dicts:
        # get 'Abiotrophia_defectiva_ATCC_49176' from 'Abiotrophia_defectiva_ATCC_49176.xml'
        folder_name = dict['name'].split('.')[0]
        # create a sub dir for every species
        temp_output_dir = os.path.join(local_output_dir,folder_name)
        # check if exists 
        if not os.path.exists(temp_output_dir):
            os.mkdir(temp_output_dir)
        # dict['name'] will be in .xml format
        file_path = os.path.join(temp_output_dir, dict['name'])
        with open(file_path, 'w') as f:
            url = dict['url']
            r = requests.get(f'{url}?raw=true')
            f.write(r.text)
        logging.info('%s written successfully.', folder_name)


def port_AGORA():
    """port different species in AGORA model

    """

    # load configure dict
    model_source = Sources.AGORA
    cur_info = basic_info_dict[Sources.AGORA]

    # logging.info(f""""
    #             ##################################
    #             #### start downloading AGORA  ####
    #             ##################################  
    # """)

    # git_download_from_direcotry(cur_info['dir_url'], cur_info['local_output_dir'])

    # logging.info(f""""
    #             ##################################
    #             #### finish downloading AGORA ####
    #             ##################################  
    # """)

    logging.info(f""""
                ##################################
                ####   start porting AGORA    ####
                ##################################  
    """)

    # load the xml files
    forder_paths = [os.path.join(cur_info['local_output_dir'], x) for x in os.listdir(cur_info['local_output_dir'])]
    print(forder_paths)

    for forder_path in forder_paths:
        file_path = os.path.join(forder_path, "{}.xml".format(forder_path.split('/')[-1]))
        
        # may encounter "'./testdata/AGORA/.DS_Store/.DS_Store.xml' does not exist"
        if not os.path.exists(file_path):
            continue
        
        # read xml by cobra
        model = cobra.io.read_sbml_model(file_path)

        ##################################
        #####     port metabolites   #####
        ##################################

        myCpds = [port_metabolite(metabolite, model_source) for metabolite in model.metabolites]

        logging.info(f'Before decompartmentalization, there are {len(myCpds)} compounds')
        # remove duplicated compounds
        myCpds = remove_duplicate_cpd(myCpds)

        logging.info(f'After decompartmentalization, there are {len(myCpds)} compounds left')

        myCpds = fetch_AGORA_GEM_identifiers(myCpds,json_path = './jms/data/staged/vmh.json',overwrite = True)


        ##################################
        #####     port reactions     #####
        ##################################
        myRxns = [port_reaction(reaction, model_source)for reaction in model.reactions]
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
        myPathways = [port_pathway(pathway, model_source, None)for pathway in model.groups]
        logging.info(
            f'There are {len(myPathways)} groups or pathways in the model')
        # retain the valid reactions in list of pathway
        myPathways = retain_valid_Rxns_in_Pathways(myPathways, myRxns)

        ##################################
        #####     build model        #####
        ##################################
        # metabolicModel to export
        MM = MetabolicModel()
        MM.meta_data = cur_info['metadata']
        species = forder_path.split('/')[-1]
        MM.id = f'{species}' 
        MM.list_of_pathways = [P.serialize() for P in myPathways]
        MM.list_of_reactions = [R.serialize() for R in myRxns]
        MM.list_of_compounds = [C.serialize() for C in myCpds]

        #################################
        ####      output model      #####
        #################################
        # Write pickle file
        export_pickle(os.path.join(forder_path, f'{MM.id}.pickle'), MM)
        logging.info(f'Export pickle file')

        # Write json file
        export_json(os.path.join(forder_path, f'{MM.id}.json'), MM)
        logging.info(f'Export json file')

        # Write dataframe
        export_table(os.path.join(
            forder_path, f'{MM.id}_list_of_compounds.csv'), MM, 'list_of_compounds')
        export_table(os.path.join(
            forder_path, f'{MM.id}_list_of_reactions.csv'), MM, 'list_of_reactions')
        export_table(os.path.join(
            forder_path, f'{MM.id}_list_of_pathways.csv'), MM, 'list_of_pathways')

        logging.info(f'Export a table of the list of compounds')

        logging.info(f""""
                    ##################################
                    #### finish porting {species} ####
                    ##################################  
        """)
