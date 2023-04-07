import os
import cobra
import logging

from jms.port.port_config import *
from jms.port.port_utils import *

from metDataModel.core import MetabolicModel

# Logger Configuration
formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logging.basicConfig(level=logging.INFO, handlers=[handler])

def port_GAPSEQ():
    """port different species in GAPSEQ model

    """

    # load configure dict
    model_source = Sources.GAPSEQ
    cur_info = basic_info_dict[model_source]

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

        model_name = "{}-gapseq".format(forder_path.split('/')[-1])

        myCpds = [port_metabolite(metabolite, model_source, model_name) for metabolite in model.metabolites]

        logging.info(f'Before decompartmentalization, there are {len(myCpds)} compounds')
        # remove duplicated compounds
        myCpds = remove_duplicate_cpd(myCpds)

        logging.info(f'After decompartmentalization, there are {len(myCpds)} compounds left')

        ##################################
        #####     port reactions     #####
        ##################################
        myRxns = [port_reaction(reaction, model_source)for reaction in model.reactions if reaction]
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
        MM.meta_data = cur_info['meta_data']
        MM.meta_data['species'] = model_name.split('-')[0]
        MM.id = f'az_{model_name}_{today}' 
        MM.list_of_pathways = [P.serialize() for P in myPathways]
        MM.list_of_reactions = [R.__dict__  for R in myRxns]
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
        species = forder_path.split('/')[-1]
        logging.info(f""""
                    ##################################
                    #### finish porting {species} ####
                    ##################################  
        """)