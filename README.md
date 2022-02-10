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
# A LC-MS peak
{
    'id_number': 555,
    'mz': 133.0970, 
    'apex': 654, 
    'height': 14388.0, 
    'left_base': 648, 
    'right_base': 655, 
}

# An empirical compound:
{
    "neutral_formula_mass": 268.08077, 
    "neutral_formula": C10H12N4O5,
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

The App Engine web tool is in a different repository.

How to use
==========

to come.


Related libraries and tools
===========================

mummichog(3)

Azimuth DB: the chemical database for biology, including metabolic models

metDataModel: data models for metabolomics, used by mummichog and Azimuth DB

mass2chem: common utilities in interpreting mass spectrometry data, annotation

Asari: trackable and scalable LC-MS data preprocessing


Citation
========
To come.

