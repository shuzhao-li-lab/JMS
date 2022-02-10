'''
# parse HMDB ver 4, single XML files for all metabolites
# Retrieved from https://hmdb.ca/downloads
# Using Python 3 and lxml
# lxml is not in current Python 3 package. 
# It has similar API as ElementTree but faster. To install:
# pip install lxml

In [1]: # cut HMDB metabolite file into 5 pieces                                               

In [2]: infile = 'xml/hmdb_metabolites.xml'                                                    

In [3]: w = open(infile).read()                                                                

In [4]: len(w)                                                                                 
Out[4]: 4077405956

In [5]: len(w)/5                                                                               
Out[5]: 815481191.2

In [6]: w[-30:]                                                                                
Out[6]: 'ations>\n</metabolite>\n</hmdb>\n'

In [7]: w[:30]                                                                                 
Out[7]: '<?xml version="1.0" encoding="'

In [8]: w[:300]                                                                                
Out[8]: '<?xml version="1.0" encoding="UTF-8"?>\n<hmdb xmlns="http://www.hmdb.ca">\n<metabolite>\n  <version>4.0</version>\n  <creation_date>2005-11-16 15:48:42 UTC</creation_date>\n  <update_date>2020-04-23 20:53:36 UTC</update_date>\n  <accession>HMDB0000001</accession>\n  <status>quantified</status>\n  <secondary'

In [9]: w[-300:]                                                                               
Out[9]: 'ctrometry for metabolomic studies of human urine. Anal Chem. 2012 Feb 21;84(4):1994-2001. doi: 10.1021/ac2030738. Epub 2012 Feb 1.</reference_text>\n      <pubmed_id>22409530</pubmed_id>\n    </reference>\n  </general_references>\n  <protein_associations>\n  </protein_associations>\n</metabolite>\n</hmdb>\n'

In [10]: block = 815481191                                                                     

In [11]:                                                                                       

In [11]: w[:300].find('<metabolite>')                                                          
Out[11]: 73

In [12]: w[:300][73]                                                                           
Out[12]: '<'

In [13]: w[:300][72]                                                                           
Out[13]: '\n'

In [14]: w[block:].find('<metabolite>')                                                        
Out[14]: 14578

In [15]: part1 = w[: block+14578]                                                              

In [16]: part1[-11:]                                                                           
Out[16]: 'etabolite>\n'

In [17]: part1[-100:]                                                                          
Out[17]: 'eferences>\n  </general_references>\n  <protein_associations>\n  </protein_associations>\n</metabolite>\n'

In [18]: tail = '</hmdb>\n'                                                                    

In [19]: with open('HMDB4_part1.xml', 'w') as O: 
    ...:     O.write(part1 + tail) 
    ...:                                                                                       

In [20]: del part1                                                                             

In [21]: w = w[block+14578:]                                                                   

In [22]:                                                                                       

In [22]: w[:80]                                                                                
Out[22]: '<metabolite>\n  <version>4.0</version>\n  <creation_date>2012-09-13 11:52:34 UTC</'

In [23]: w[block:].find('<metabolite>')                                                        
Out[23]: 7241

In [24]: part2 = w[: block+7241]                                                               

In [25]: w = w[block+7241:]                                                                    

In [26]: with open('HMDB4_part2.xml', 'w') as O: 
    ...:     O.write(part2+tail) 
    ...:                                                                                       

In [27]: del part2                                                                             

In [28]:                                                                                       

In [28]: w[block:].find('<metabolite>')                                                        
Out[28]: 4790

In [29]: head = '<?xml version="1.0" encoding="UTF-8"?>\n<hmdb xmlns="http://www.hmdb.ca">\n'  

In [30]: part3 = w[:block+4790]                                                                

In [31]: w = w[block+4790:]                                                                    

In [32]: with open('HMDB4_part3.xml', 'w') as O: 
    ...:     O.write( head + part3 + tail ) 
    ...:                                                                                       

In [33]: del part3                                                                             

In [34]:                                                                                       

In [34]: w[block:].find('<metabolite>')                                                        
Out[34]: 11469

In [35]: part4 = w[:block+11469]                                                               

In [36]: w = w[block+11469:]                                                                   

In [37]:                                                                                       

In [37]: with open('HMDB4_part4.xml', 'w') as O: 
    ...:     O.write( head + part4 + tail ) 
    ...:                                                                                       

In [38]: del part4                                                                             

In [39]:                                                                                       

In [39]: with open('HMDB4_part5.xml', 'w') as O: 
    ...:     O.write( head + w ) 
    ...:                                                                                       

In [40]: del w                                                                                 

In [41]:                                                                                       

In [41]: # fix part2                                                                           

In [42]: part2 = open('HMDB4_part2.xml').read()                                                

In [43]: with open('HMDB4_part2.xml', 'w') as O: 
    ...:     O.write( head + part2 ) 
    ...:                                                                                       

In [44]: del part2                                                                             

In [45]:     

# parsing using this script

# This took ~5 seconds on my i5 laptop, infile is 1.2 GB.
# Larger files will have to be partitioned, depending on memeory limit.
In [1]: from parse_hmdb4 import *                                                              
In [2]: myList = parse_a_file(infile, wanted)                                                  
Found 25411 entries.
In [3]: myList[0]                                                                              
Out[3]: 
{'accession': 'HMDB0000001',
 'name': '1-Methylhistidine',
 'chemical_formula': 'C7H11N3O2',
 'monisotopic_moleculate_weight': '',
 'iupac_name': '(2S)-2-amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid',
 'traditional_iupac': '1 methylhistidine',
 'cas_registry_number': '332-80-9',
 'smiles': 'CN1C=NC(C[C@H](N)C(O)=O)=C1',
 'inchi': 'InChI=1S/C7H11N3O2/c1-10-3-5(9-4-10)2-6(8)7(11)12/h3-4,6H,2,8H2,1H3,(H,11,12)/t6-/m0/s1',
 'inchikey': 'BRMWTNUJHUMWMS-LURJTMIESA-N',
 'pathways': '',
 'normal_concentrations': '',
 'abnormal_concentrations': '',
 'diseases': '',
 'drugbank_id': 'DB04151',
 'drugbank_metabolite_id': '',
 'phenol_explorer_compound_id': '',
 'phenol_explorer_metabolite_id': '',
 'foodb_id': 'FDB093588',
 'knapsack_id': '',
 'chemspider_id': '83153',
 'kegg_id': 'C01152',
 'biocyc_id': '',
 'bigg_id': '',
 'wikipidia': '',
 'nugowiki': '',
 'metagene': '',
 'metlin_id': '3741',
 'pubchem_compound_id': '92105',
 'het_id': '',
 'chebi_id': '50599',
 'protein_associations': ''}
In [40]: with open('HMDB0000031.xml', 'w') as O: 
    ...:     O.write(str(ET.tostring(x), 'utf-8')) 
    ...:                                                                                       
# Alternatively, cut the big file into 5 parts, then parse:
In [2]: myResults = []                                                                         
In [3]: wanted = ['accession', 'name', 'chemical_formula', 'monisotopic_molecular_weight', 'keg
   ...: g_id', 'pubchem_compound_id', ]                                                        
In [4]: import os                                                                              
In [5]: files = os.listdir('parts/')                                                           
In [6]: files                                                                                  
Out[6]: 
['README.txt',
 'HMDB4_part5.xml',
 'HMDB4_part3.xml',
 'HMDB4_part2.xml',
 'HMDB4_part1.xml',
 'HMDB4_part4.xml']
In [7]: for f in files[1:]: 
   ...:     myResults += parse_a_file('parts/' + f, wanted) 
   ...:                                                                                        
Found 30770 entries.
Found 21801 entries.
Found 17478 entries.
Found 22506 entries.
Found 21667 entries.
In [8]: write_tsv(myResults, wanted, "HMDB4_compounds.tsv")                
'''

from lxml import etree as ET 

infile = 'serum_metabolites.xml'

# fields to retrieve
wanted = ['accession', 'name', 'chemical_formula', 'monisotopic_molecular_weight', 'iupac_name', 
    'traditional_iupac', 'cas_registry_number', 'smiles', 'inchi', 'inchikey', 'pathways', 
    'normal_concentrations', 'abnormal_concentrations', 'diseases', 'drugbank_id', 
    'drugbank_metabolite_id', 'phenol_explorer_compound_id', 'phenol_explorer_metabolite_id', 
    'foodb_id', 'knapsack_id', 'chemspider_id', 'kegg_id', 'biocyc_id', 'bigg_id', 'wikipidia', 
    'nugowiki', 'metagene', 'metlin_id', 'pubchem_compound_id', 'het_id', 'chebi_id', 'protein_associations']

def extract(obj_metabolite, x):
    try:
        y = obj_metabolite.find(x).text.strip()
        if y:
            return y
        else: return ''
    except AttributeError: 
        return ''

def extract_dict(obj_metabolite, wanted, prefix):
    result = {}
    for x in wanted:
        result[x] = extract(obj_metabolite, prefix+x)
    return result

def parse_a_file(f, wanted, prefix='{http://www.hmdb.ca}'):
    '''
    Input
    =====
    f: HMDB 4 XML file
    wanted: fields to retrieve
    prefix: prefix that is used in the XML file
    Return
    ======
    return as a list of dictionaries. 
    '''
    tree = ET.parse(f)
    root = tree.getroot()
    print("Found %d entries." %len(root))

    results = []
    for child in root:
        results.append(extract_dict(child, wanted, prefix=prefix))

    return results

def write_tsv(results, wanted, outfile):
    s = '\t'.join(wanted) + '\n'
    for R in results:
        s += '\t'.join([R[x] for x in wanted]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)

if __name__ == '__main__':
    # wanted = ['accession', 'name', 'chemical_formula', 'monisotopic_molecular_weight']
    infile = 'serum_metabolites.xml'
    # infile= 'urine_metabolites.xml'
    # infile= 'feces_metabolites.xml'

    write_tsv(
        parse_a_file(infile, wanted),
        wanted, 
        "parsed_"+infile.replace(".xml", ".tsv")
    )
