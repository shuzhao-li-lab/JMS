"""
Function used to port different genome-scale metabolite models into metDataModel.

Yuanye Chi 2023-03-22
"""

from jms.port.port_config import *
from jms.port.port_CUT import port_CUT
from jms.port.port_AGORA import port_AGORA
    

def port(model_source):

    if model_source == Sources.CUT:
        for species in basic_info_dict[Sources.CUT].keys():
            port_CUT(species)
    elif model_source == Sources.AGORA:
        port_AGORA()
        


 