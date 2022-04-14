JMS Documenetation
==================

dbStructures contain two main classes: ExperimentalEcpdDatabase, knownCompoundDatabase

Example compound in KCD.mass_indexed_compounds:
```
    ('C8H13N2O2_169.097154',
    {'interim_id': 'C8H13N2O2_169.097154',
    'neutral_formula': 'C8H13N2O2',
    'neutral_formula_mass': 169.097154087,
    'compounds': [{'primary_id': 'HMDB0062696',
        'primary_db': 'HMDB',
        'name': 'Pyridoxaminium(1+)',
        'neutral_formula': 'C8H13N2O2',
        'neutral_formula_mass': 169.097154087,
        'SMILES': 'CC1=C(O)C(C[NH3+])=C(CO)C=N1',
        'inchikey': 'NHZMQXZHNVQTQA-UHFFFAOYSA-O',
        'other_ids': {'PubChem': '25245492', 'KEGG': '', 'ChEBI': '57761'}}]})
```

Example search result:
```
        peak_result_dict can have one empCpd matched to multiple DB entries, e.g.
            {'peak': {'parent_masstrack_id': 240,
            'mz': 119.12584100277608,
            'id_number': 'F201'},
            'interim_id': 20,
            'epd_ion_relation': '13C/12C',
            'list_matches': [('C6H15NO_117.115364', 'M+H[1+]', 2),
            ('C6H13N_99.104799', 'M+H2O+H[1+]', 1)]}
```

Conversion from GEM to JSON
===========================






