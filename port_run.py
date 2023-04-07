"""
Function used to port different genome-scale metabolite models into metDataModel.

Yuanye Chi 2023-03-22
"""
import argparse
import logging

from jms.port.port_config import Sources, basic_info_dict
from jms.port.port_CUT import port_CUT
from jms.port.port_AGORA import port_AGORA
from jms.port.port_GAPSEQ import port_GAPSEQ


if __name__ == '__main__':
    # jm.port(Sources.CUT)
    # jm.port(Sources.AGORA)
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", help="choose a mdoel you want to port in.", choices=["CUT", "AGORA","GAPSEQ"], required=True)
    args = parser.parse_args()

    logging.info(args)

    if args.model == "CUT":
        for species in basic_info_dict[Sources.CUT].keys():
            port_CUT(species)
    elif args.model == "AGORA":
        port_AGORA()
    elif args.model == "GAPSEQ":
        port_GAPSEQ()