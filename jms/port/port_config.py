from enum import Enum

from datetime import datetime
today = str(datetime.today()).split(" ")[0]


class Sources(Enum):
    CUT = 1
    AGORA = 2


basic_info_dict = {
    Sources.CUT: {
        'human': {
            # The path you intended to save
            'output_dir': './model_output/HumanGEM/',
            'name_xml': 'Human-GEM.xml',
            'github_xml_path': 'https://github.com/SysBioChalmers/Human-GEM/blob/main/model/Human-GEM.xml',
            'note': 'Human-GEM compartmentalized, with genes and ECs.',
            'pathway_source': 'Human-GEM v1.14.0',
            'metadata': {
                'species': 'human',
                'version': '',
                'sources': [f'https://github.com/SysBioChalmers/Human-GEM, retrieved {today}'],
                'status': '',
                'last_update': today,  #
                'note': 'Human-GEM compartmentalized, with genes and ECs.',
            }
        },
        'mouse': {
            # The path you intended to save
            'output_dir': './model_output/MouseGEM/',
            'name_xml': 'Mouse-GEM.xml',
            'github_xml_path': 'https://github.com/SysBioChalmers/Mouse-GEM/blob/main/model/Mouse-GEM.xml',
            'note': 'Mouse-GEM compartmentalized, with genes and ECs.',
            'pathway_source': 'Mouse-GEM v1.4.0',
            'metadata': {
                'species': 'mouse',
                'version': '',
                'sources': [f'https://github.com/SysBioChalmers/Mouse-GEM, retrieved {today}'],
                'status': '',
                'last_update': today,  #
                'note': 'Mouse-GEM compartmentalized, with genes and ECs.',
            }
        },
        'rat': {
            # The path you intended to save
            'output_dir': './model_output/RatGEM/',
            'name_xml': 'Rat-GEM.xml',
            'github_xml_path': 'https://github.com/SysBioChalmers/Rat-GEM/blob/main/model/Rat-GEM.xml',
            'note': 'Rat-GEM compartmentalized, with genes and ECs.',
            'pathway_source': 'Rat-GEM v1.4.0',
            'metadata': {
                'species': 'rat',
                'version': '',
                'sources': [f'https://github.com/SysBioChalmers/Rat-GEM, retrieved {today}'],
                'status': '',
                'last_update': today,  #
                'note': 'Rat-GEM compartmentalized, with genes and ECs.',
            }
        },
        'worm': {
            # The path you intended to save
            'output_dir': './model_output/WormGEM/',
            'name_xml': 'Worm-GEM.xml',
            'github_xml_path': 'https://github.com/SysBioChalmers/Worm-GEM/blob/main/model/Worm-GEM.xml',
            'note': 'Worm-GEM compartmentalized, with genes and ECs.',
            'pathway_source': 'Worm-GEM v1.4.0',
            'metadata': {
                'species': 'worm',
                'version': '',
                'sources': [f'https://github.com/SysBioChalmers/Worm-GEM, retrieved {today}'],
                'status': '',
                'last_update': today,  #
                'note': 'Worm-GEM compartmentalized, with genes and ECs.',
            }
        },
        'yeast': {
            # The path you intended to save
            'output_dir': './model_output/YeastGEM/',
            'name_xml': 'Yeast-GEM.xml',
            'github_xml_path': 'https://github.com/SysBioChalmers/yeast-GEM/blob/main/model/yeast-GEM.xml',
            'note': 'Yeast-GEM compartmentalized, with genes and ECs.',
            'pathway_source': 'Yeast-GEM v8.6.2',
            'metadata': {
                'species': 'yeast',
                'version': '',
                'sources': [f'https://github.com/SysBioChalmers/yeast-GEM, retrieved {today}'],
                'status': '',
                'last_update': today,  #
                'note': 'Yeast-GEM compartmentalized, with genes and ECs.',
            }
        },
        'zebrafish': {
            # The path you intended to save
            'output_dir': './model_output/ZebrafishGEM/',
            'name_xml': 'Zebrafish-GEM.xml',
            'github_xml_path': 'https://github.com/SysBioChalmers/Zebrafish-GEM/blob/main/model/Zebrafish-GEM.xml',
            'note': 'Zebrafish-GEM compartmentalized, with genes and ECs.',
            'pathway_source': 'Zebrafish-GEM v1.3.0',
            'metadata': {
                'species': 'zebrafish',
                'version': '',
                'sources': [f'https://github.com/SysBioChalmers/Zebrafish-GEM, retrieved {today}'],
                'status': '',
                'last_update': today,  #
                'note': 'Zebrafish-GEM compartmentalized, with genes and ECs.',
            }
        },
        'fruitfly': {
            # The path you intended to save
            'output_dir': './model_output/ZebrafishGEM/',
            'name_xml': 'Fruitfly-GEM.xml',
            'github_xml_path': 'https://github.com/SysBioChalmers/Fruitfly-GEM/blob/main/model/Fruitfly-GEM.xml',
            'note': 'Fruitfly-GEM compartmentalized, with genes and ECs.',
            'pathway_source': 'Fruitfly-GEM v7.7.7',
            'metadata': {
                'species': 'fruitfly',
                'version': '',
                'sources': [f'https://github.com/SysBioChalmers/Fruitfly-GEM, retrieved {today}'],
                'status': '',
                'last_update': today,  #
                'note': 'Fruitfly-GEM compartmentalized, with genes and ECs.',
            }
        }
    },
    Sources.AGORA: {
        # url to get the AGORA 
        'dir_url': "https://api.github.com/repos/VirtualMetabolicHuman/AGORA/contents/CurrentVersion/AGORA_1_03/AGORA_1_03_With_Mucins_sbml",   
        'local_output_dir': "./testdata/AGORA",
        'metadata': {
                    'species': '',
                    'version': 'AGORA_1_03_With_Mucins_sbml',
                    'sources': [f'AGORA, retrieved {today}'], #
                    'status': '',
                    'last_update': today,  #
                    'note': f'AGORA cloned from https://github.com/VirtualMetabolicHuman, retrieved from {today}\ .'
                }
    }
}
