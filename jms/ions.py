'''
Functions related to isotopes and adductions.
Using mass2chem.

compute_adducts_formulae(304.2038, 'C19H28O3',  mode='pos', primary_only=True)
Out[5]:
[[304.203251, 'M[1+]', 'C19H28O3'],
 [305.21107646677, 'M+H[1+]', 'C19H29O3'],
 [327.19307646677, 'M+Na[1+]', 'C19H28NaO3'],
 [323.22167646677, 'M+H2O+H[1+]', 'C19H31O4']]

'''


from mass2chem.formula import compute_adducts_formulae



