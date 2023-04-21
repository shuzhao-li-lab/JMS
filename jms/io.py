'''
Input/Output functions

>>> import sys
>>> sys.path.append('jms')
>>> from jms.io import read_hmdb_to_empCpds
>>> f = '/Users/shuzhao/li.github/mass2chem/mass2chem/source_data/HMDB4_compounds.tsv'
>>> list_empCpds = read_hmdb_to_empCpds(f)
>>> list_empCpds[5559]
{'interim_id': 'C13H18N4O6_326.122634', 'neutral_formula': 'C13H18N4O6', 'neutral_formula_mass': 326.122634328, 
'compounds': [{'id': 'HMDB0003826', 'name': '6,7-Dimethyl-8-(1-D-ribityl)lumazine'}]}
>>> import json
>>> with open('jms/data/compounds/hmdb4_compounds.json', 'w', encoding='utf-8') as f:
...   json.dump(list_empCpds, f, ensure_ascii=False, indent=2)
... 
>>> len(list_empCpds)
11937

# example of a Compound:
        {'primary_id': 'HMDB0042472', 'primary_db': 'HMDB', 
        'name': 'TG(14:0/22:1(13Z)/22:4(7Z,10Z,13Z,16Z))', 
        'neutral_formula': 'C61H108O6', 'neutral_formula_mass': 936.814591198, 
        'SMILES': '[H][C@](COC(=O)CCCCCCCCCCCCC)(COC(=O)CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC)OC(=O)CCCCCCCCCCC\\C=C/CCCCCCCC', 
        'inchikey': 'NDYSJNXSSRKHCQ-LHDKJNFNSA-N', 
        'other_ids': {'PubChem': '131753609', 'KEGG': '', 'ChEBI': ''}}
'''

peak_attribute_dict = {
    'm/z': 'mz',
    'mz': 'mz',
    'apex': 'rtime',
    # ...,
}

def read_table_to_peaks(infile, 
                        has_header=True, mz_col=1, rtime_col=2, intensity=(11,21), feature_id=None,
                        full_extract=True, max_col=21,
                        delimiter='\t'):
    '''
    Read a text feature table, and 
    return list of peaks, e.g. [ {
        'id_number': 555,
        'mz': 133.0970, 
        'apex': 654, 
        'height': 14388.0, 
        }, ... ]
    
    full_extract: to keep all fields in output as strings, only if has_header.
    feature_id: if None, create id for each peak/feature.
    '''
    def _make_id(ii, mz, rt):
        return 'F' + str(ii) + '_' + str(round(mz, 6)) + '@' + str(round(rt, 2))

    list_peaks = []
    w = open(infile).readlines()
    if has_header:
        header = w[0].rstrip().split(delimiter)[:max_col]
        w = w[1:]
    ii = 0
    for line in w:
        a = line.split(delimiter)[:max_col]   # not rstrip, so trailing EOL will be carried forward
        mz, rt = float(a[mz_col]), float(a[rtime_col])
        intensities = [float(x) for x in a[intensity[0]: intensity[1]]]
        if feature_id != None:
            fid = a[feature_id].strip()
        else:
            ii += 1
            fid = _make_id(ii, mz, rt)
        peak = {'id_number': fid, 'mz': mz, 'rtime': rt, 'apex': rt, 'representative_intensity': sum(intensities) / len(intensities)}
        if has_header and full_extract:
            # will remove redundant fields
            peak2 = dict(zip(header, a))
            peak2.update(peak)
            peak = peak2

        list_peaks.append( peak )

    return list_peaks


def read_tsv_hmdb_to_empCpds(infile, delimiter='\t'):
    '''
    Previously parsed HMDB4 in short tsv format.
    Ref:
    https://github.com/shuzhao-li/JMS/wiki/Scripting-to-create-a-dictionary-from-HMDB-data
    https://github.com/shuzhao-li/JMS/blob/main/jms/utils/hmdb.py
    A copy of HMDB4_compounds.tsv is under
    https://github.com/shuzhao-li/mass2chem/tree/master/mass2chem/lib

    In [26]: hmdb[0].split('\t')                                                                                              
    Out[26]: 
    ['accession',
    'name',
    'chemical_formula',
    'monisotopic_molecular_weight',
    'kegg_id',
    'pubchem_compound_id\n']

    This is not producing exact format as in dbStructures, which should be followed when needed more elaborate.
    '''
    hmdb = open(infile).readlines()
    hdict = {}                                                                                                       
    for line in hmdb[1:]: 
        a = line.split(delimiter) 
        if a[3].strip(): 
            hdict[a[0]] = a 
    h0 = [x.split('\t')[0] for x in hmdb[1:]] 
    massDict_hmdb = {}
    for x in h0:  
        try: 
            a = hdict[x]  
            k = a[2]+ '_' + str(round(float(a[3]),6))  # ensuring identical formula and mass
            if k in massDict_hmdb:  
                massDict_hmdb[k].append( a )  
            else:  
                massDict_hmdb[k] = [a]  
        except KeyError: 
            print(x) 
    list_empCpds = []
    for k,v in massDict_hmdb.items():
        compounds = []
        for cpd in v:
            compounds.append({
                "id": cpd[0], "name": cpd[1]
            })
        list_empCpds.append(
            {
                "interim_id": k,
                "neutral_formula": v[0][2],
                "neutral_formula_mass": float(v[0][3]),
                "compounds": compounds,
            }
        )
    return list_empCpds


