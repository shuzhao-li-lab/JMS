"""
Function used to port different genome-scale metabolite models into metDataModel.

Yuanye Chi 2023-03-22
"""
import sys
import os
from datetime import datetime

import cobra # https://cobrapy.readthedocs.io/en/latest/io.html#SBML

import requests

from mass2chem.formula import *

from jms.formula import *
from jms.port.port_utils import *
from jms.port.port_config import *
from jms.port.port_CUT import port_CUT
from jms.port.port_AGORA import port_AGORA

from metDataModel.core import Compound, Reaction, Pathway, MetabolicModel
    

def port(model_source):

    if model_source == Sources.CUT:
        for species in basic_info_dict[Sources.CUT].keys():
            port_CUT(species)
    elif model_source == Sources.AGORA:
        port_AGORA()
        


 


if __name__ == '__main__':
    port(Sources.CUT)