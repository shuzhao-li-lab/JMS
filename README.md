Json's Metabolite Services (JMS)
================================

A Python library for 
- mapping identifiers between genome scale metabolic models and metabolite databases
- reusable data structures for peaks, compounds and indexed data stores
- efficient mass and empirical compound search functions

As the name suggests, JSON or similar Python dictionaries are used widely in this library.
The definition of these data structures was a result of patterning our scientific work, 
and we encourage reuse and community contributions.
A few examples:

```
// A LC-MS peak:
{
    'id_number': 555,
    'mz': 133.0970, 
    'rtime': 654, 
    'height': 14388.0, 
    'left_base': 648, 
    'right_base': 662, 
}

// A compound (i.e. metabolite):
{
    'primary_id': HMDB0000195,
    'primary_db': 'HMDB',
    'name': 'Inosine',
    "neutral_formula": C10H12N4O5,
    "neutral_formula_mass": 268.08077, 
    'SMILES': 'OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=NC2=C(O)N=CN=C12', 
    'inchikey': 'UGQMRVRMYYASKQ-KQYNXXCUSA-N',
    'LogP': -2.10,
    'other_ids': {'PubChem': '6021',
                  'KEGG': 'C00294',
                  'ChEBI': '17596',
                  'MetaNetX': 'MNXM1103335',
                  },
}

// An empirical compound:
{
    "neutral_formula_mass": 268.08077, 
    "neutral_formula": C10H12N4O5,
    "biological_ion": '',
    "biological_ion_charge": None,
    "alternative_formulas": [],
    "interim_id": C10H12N4O5_268.08077,
    "identity": [
            {'compounds': ['HMDB0000195'], 'names': ['Inosine'], 'score': 0.6, 'probability': null},
            {'compounds': ['HMDB0000195', 'HMDB0000481'], 'names': ['Inosine', 'Allopurinol riboside'], 'score': 0.1, 'probability': null},
            {'compounds': ['HMDB0000481'], 'names': ['Allopurinol riboside'], 'score': 0.1, 'probability': null},
            {'compounds': ['HMDB0003040''], 'names': ['Arabinosylhypoxanthine'], 'score': 0.05, 'probability': null},
            ],
    "MS1_pseudo_Spectra": [
            {'feature_id': 'FT1705', 'mz': 269.0878, 'rtime': 99.90, 'charged_formula': '', 'ion_relation': 'M+H[1+]'},
            {'feature_id': 'FT1876', 'mz': 291.0697, 'rtime': 99.53, 'charged_formula': '', 'ion_relation': 'M+Na[1+]'},
            {'feature_id': 'FT1721', 'mz': 270.0912, 'rtime': 99.91, 'charged_formula': '', 'ion_relation': 'M(C13)+H[1+]'},
            {'feature_id': 'FT1993', 'mz': 307.0436, 'rtime': 99.79, 'charged_formula': '', 'ion_relation': 'M+K[1+]'},
            ],
    "MS2_Spectra": ['AZ0000711', 'AZ0002101'],
    "Database_referred": ["Azimuth", "HMDB", "MONA"],
}
```

These examples show a core set of attributes. 
Users can extend to include other attributes that meet their own project needs.
With consistent core concepts, some simple mapping dictionaries will produce interoperability between projects. 

In empirical compound, MS1_pseudo_Spectra are for adduction ions (A-ions); the biological_ion is usually different (B-ions).

Packaging
=========

The PyPi package excludes subdirectories in /data.


The App Engine web tool is in a different repository.

How to use
==========

See notebooks under docs/ for demonstration.

Example code of annotating a feature table against HMDB:

```
import json
from jms.dbStructures import annotate_peaks_against_kcds
from jms.io import read_table_to_peaks

list_compounds = json.load(open('jms/data/compounds/list_compounds_HMDB4.json'))
mydata = read_table_to_peaks('pos_cmap_feature_table.csv', ',')

annotate_peaks_against_kcds(mydata, list_compounds, 
                                export_file_name_prefix='jms_annotated_',
                                mode='pos',  mz_tolerance_ppm=5)
```

Related libraries and tools
===========================

metDataModel: data models for metabolomics

mass2chem: common utilities in interpreting mass spectrometry data, annotation

Asari: trackable and scalable LC-MS data preprocessing

mummichog(3)

Azimuth DB: the chemical database for biology, including metabolic models


Citation
========
To come.

