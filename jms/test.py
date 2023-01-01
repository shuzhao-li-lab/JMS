'''
# decompress the HMDB4 data if needed
➜  JMS git:(main) ✗ xz -d jms/data/compounds/list_compounds_HMDB4.json.xz 
# run test from top dir
➜  JMS git:(main) ✗ python3 -m jms.test 
'''

import json
from jms.dbStructures import annotate_peaks_against_kcds
from jms.io import read_table_to_peaks

list_compounds = json.load(open('jms/data/compounds/list_compounds_HMDB4.json'))
mydata = read_table_to_peaks('testdata/full_Feature_table.tsv', '\t')

annotate_peaks_against_kcds(mydata, list_compounds, 
                                export_file_name_prefix='jms_annotated_',
                                mode='pos',  mz_tolerance_ppm=5)

