'''
Use lists from mass2chem.

It's currently combined from:

Keller, B.O.; Sui, J.; Young, A.B.; Whittal, R.M. 
Interferences and contaminants encountered in modern mass spectrometry.
Analytica Chimica Acta (Review/tutorial, Special Issue on Mass Spectrometry), 2008.

Ralf J. M. Weber, Eva Li, Jonathan Bruty, Shan He & Mark R. Viant (2012). 
MaConDa: a publicly accessible Mass spectrometry Contaminants Database. 
Bioinformatics, doi:10.1093/bioinformatics/bts527.

# Example of wrangling:

list_contaminants_pos, list_contaminants_neg = [], []
dict_contaminants_pos, dict_contaminants_neg = {}, {}

for x in contaminants_pos:
    list_contaminants_pos.append({
        'id': x[1],
        'mz': x[0], # - PROTON,
        'rtime': 0,
        # 'full': x,
    })
    dict_contaminants_pos[x[1]] = x
    
for x in contaminants_neg:
    list_contaminants_neg.append({
        'id': x[1],
        'mz': x[0], # + PROTON,
        'rtime': 0,
        # 'full': x,
    })  
    dict_contaminants_neg[x[1]] = x

print( len(list_contaminants_pos), len(list_contaminants_neg),
      list_contaminants_neg[3] )
      

>>> 
1297 39 {'id': '65Cu(I)_[M+I]-_318.737295', 'mz': 318.7372946, 'rtime': 0}

'''

from mass2chem.lib.LCMS_contaminants import contaminants_pos, contaminants_neg
