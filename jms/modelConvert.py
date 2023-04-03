'''
Conversion of genome scale metabolic models to JSON formats,and to indexed Python internal objects.

Started in mummichog 2.6, but moving to this repo and will continue with mummichog 3.

'''

def convert_json_model(jmodel):
    '''
    Convert a JSON style, concise genome scale metabolic models to mummichog style with internal indexing.
    All compound identifiers and other identifiers are expected to be contained within the model.
    Returns model as a dictionary with as in mummichog style.

    >>> jmodel.keys()
        dict_keys(['id', 'list_of_reactions', 'list_of_compounds', 'list_of_pathways', 'meta_data'])
    >>> jmodel['list_of_compounds'][772]
    {'id': 'ebastineoh', 'name': 'Hydroxylated Ebastine', 'identifiers': {'inchi': 'InChI=1S/C32H39NO3/c1-32(2,24-34)28-17-15-25(16-18-28)30(35)14-9-21-33-22-19-29(20-23-33)36-31(26-10-5-3-6-11-26)27-12-7-4-8-13-27/h3-8,10-13,15-18,29,31,34H,9,14,19-24H2,1-2H3'}, 'neutral_formula': '', 'charge': 0, 'charged_formula': 'C32H39NO3', 'neutral_mono_mass': 0.0, 'SMILES': '', 'inchi': ''}
    >>> jmodel['list_of_pathways'][2]
    {'id': 'group3', 'name': 'Aminoacyl-tRNA biosynthesis', 
    'list_of_reactions': ['HMR_5130', 'HMR_5131', 'HMR_5132', 'HMR_5133', 'HMR_5134', 'HMR_5135', 
            'HMR_5136', 'HMR_5137', 'HMR_5138', 'HMR_5139', 'HMR_5140', 'HMR_5141', 'HMR_5142', 'HMR_5143', 
            'HMR_5144', 'HMR_5145', 'HMR_5146', 'HMR_5147', 'HMR_5148', 'HMR_5149', 'HMR_5150']}
    >>> 
    >>> jmodel['list_of_reactions'][22]
    {'id': '24_25DHVITD3tm', 'reactants': ['2425dhvitd3'], 'products': ['2425dhvitd3']}
    >>> mfn['metabolic_pathways'][51]
    {'cpds': ['C00024', 'C00027', 'strdnccoa', 'C00100', 'C01944', 'q10h2', 'dmhptcoa', 'etfrd', 'C03035', 'etfox', 'dmnoncoa', 'arachdcoa', 'dcsptn1coa', 'ptdcacoa', 'C00412', 'od2coa', 'hdcoa', 'C02593', 'dlnlcgcoa', 'vacccoa', 'C00010', 'C00016', 'adrncoa', 'C00154', 'lneldccoa', 'lnlccoa', 'c226coa', 'lnlncacoa', 'C00399', 'odecoa', 'hpdcacoa', 'C01352', 'eicostetcoa', 'tmndnccoa', 'arachcoa'], 
    'rxns': ['FAOXC2242046m', 'FAOXC200180m', 'FAOXC204', 'FAOXC160', 'FAOXC16080m', 'FAOXC180', 'FAOXC80', 'FAOXC140', 'FAOXC1811601m', 'FAOXC1811603m', 'FAOXC150m', 'FAOXC18480m', 'FAOXC1811602m', 'FAOXC16180m', 'FAOXC170m', 'FAOXC18280m', 'FAOXC182806m', 'FAOXC183803m', 'FAOXC2252053m', 'FAOXC2031836m', 'FAOXC11', 'FAOXC204184m', 'FAOXC226205m', 'FAOXC226', 'ETFQO', 'ETF'], 
    'ecs': ['1.3.99.3', '1.5.5.1'], 'id': 'mfn1v10path136', 'name': 'Fatty acid oxidation'}
    >>> mfn['metabolic_rxns'][2]
    {'reactants': ['C06611'], 'id': 'R05233', 'source': 'kegg_dre', 'ecs': ['1.1.1.1'], 'products': ['C06613'], 'cpds': ['C06611', 'C06613'], 'pathway': '3-Chloroacrylic acid degradation'}
    >>>    
    mfn['cpd2pathways'] = {
         'C06125': ['Glycosphingolipid metabolism', 'Glycosphingolipid biosynthesis - ganglioseries'], ...
         }
    >>> mfn['cpd_edges'][:5]
[   ['C00106', 'C00526'], ['C05977', 'CE6326'], ['C00043', 'G00174'], ['C00199', 'C00231'], ['C00014', 'C00655']]
    '''
    new = {}
    new['id'] = jmodel['id']
    new['version'] = jmodel['meta_data']['version']
    cpdDict,  dict_cpds_def = {}, {}
    for cpd in jmodel['list_of_compounds']:
        cpdDict[cpd['id']] = cpd    # convert_compound_azimuth_mummichog(cpd, hmdb_dict, kegg_dict)
        dict_cpds_def[cpd['id']] = cpd['name']
    new['Compounds'] = cpdDict
    new['dict_cpds_def'] = dict_cpds_def
    new['metabolic_rxns'] = jmodel['list_of_reactions']

    # bunch of indexing
    dict_rxns, edge2rxn, edge2enzyme = {}, {}, {}
    cpd_edges = []
    for rxn in jmodel['list_of_reactions']:
        dict_rxns[rxn['id']] = rxn
        _edges = []
        for cpd1 in rxn['reactants']:
            for cpd2 in rxn['products']:
                k = [cpd1, cpd2]
                cpd_edges.append(k)
                _edges.append(k)
                edge2rxn[','.join(sorted(k))] = rxn['id']
                enzymes = rxn.get('ecs', []) or rxn.get('enzymes', [])
                for e in enzymes:
                    edge2enzyme[','.join(sorted(k))] = e

    new['cpd_edges'] = cpd_edges
    new['edge2rxn'] = edge2rxn
    new['edge2enzyme'] = edge2enzyme
    # now do pathways
    metabolic_pathways = []
    for P in jmodel['list_of_pathways']:
        _p = {}
        _p['id'] = P['id']
        _p['name'] = P['name']
        _p['rxns'] = P['list_of_reactions']
        cpds, ecs, genes = [], [], []
        for ii in P['list_of_reactions']:
            rxn = dict_rxns[ii]
            cpds += rxn['reactants'] + rxn['products']
            enzymes = rxn.get('ecs', []) or rxn.get('enzymes', [])
            ecs += enzymes
            if 'genes' in rxn:
                genes += rxn['genes']
        _p['cpds'] = list(set(cpds))
        _p['ecs'] = list(set(ecs))
        _p['genes'] = list(set(genes))
        metabolic_pathways.append(_p)

    cpd2pathways = {}
    for _p in metabolic_pathways:
        for cpd in _p['cpds']:
            if cpd in cpd2pathways:
                cpd2pathways[cpd].append(_p['name'])
            else:
                cpd2pathways[cpd] = [_p['name']]

    new['metabolic_pathways'] = metabolic_pathways
    new['cpd2pathways'] = cpd2pathways

    return new



class metabolicNetwork:
    '''
    Simple wrapper class of a genome scale metabolic model in JSON.

    ?? keep or not ??
    '''
    def __init__(self, MetabolicModel):
        '''
        MetabolicModel : metabolic model in JSON style dictionary.
        '''
        #print_and_loginfo( "Loading metabolic network %s..." %MetabolicModel.version ) # version from metabolic model
        
        self.MetabolicModel = MetabolicModel
        self.network = self.build_network(MetabolicModel['cpd_edges'])
        
        self.version = MetabolicModel['version']
        self.Compounds = MetabolicModel['Compounds']
        self.metabolic_pathways = MetabolicModel['metabolic_pathways']
        self.dict_cpds_def = MetabolicModel['dict_cpds_def']
        self.cpd2pathways = MetabolicModel['cpd2pathways']
        self.edge2enzyme = MetabolicModel['edge2enzyme']
        self.total_cpd_list = self.network.nodes()
        
        
    def build_network(self, edges):
        return nx.from_edgelist( edges )
        

    def get_pathways(self):
        pass



class InputUserData:
    '''
    
    Changing to JSON list of features and list of epds


    '''
    
    def __init__(self, paradict, web=False):
        '''
        
        '''
        self.web = web
        self.paradict = paradict
        self.header_fields = []
        self.ListOfMassFeatures = []
        self.input_featurelist = []

        # entry point of data input
        self.read()
        self.update()
        
    def update(self):
        '''
        Update retention_time_rank and is_significant to all MassFeatures
        '''
        retention_times = [M['rtime'] for M in self.ListOfMassFeatures]
        self.max_retention_time = max(retention_times)

        self.max_mz = max([M['mz'] for M in self.ListOfMassFeatures])

        self.determine_significant_list(self.ListOfMassFeatures)
        
        


    def text_to_ListOfMassFeatures(self, textValue, delimiter='\t'):
        '''
        Column order is hard coded for now, as mz, retention_time, p_value, statistic, CompoundID_from_user

        use asari style JSON features

        '''
        def _make_id(ii, mz, rt):
            return 'F' + str(ii) + '_' + str(round(mz, 6)) + '@' + str(round(rt, 2))
        #
        lines = self.__check_redundant__( textValue.splitlines() )
        self.header_fields = lines[0].rstrip().split(delimiter)

        excluded_list = []
        for ii in range(len( lines )-1):
            y = lines[ii+1].split('\t')
            
            fid_from_user = ''
            if len(y) > 4: fid_from_user = y[4].strip()

            [mz, rtime, p_value, statistic] = [float(x) for x in y[:4]]
            
            # row_number, mz, retention_time, p_value, statistic, CompoundID_from_user
            if MASS_RANGE[0] < mz < MASS_RANGE[1]:
                # row # human-friendly, numbering from 1
                fid = _make_id(ii+1, mz, rtime)
                peak = {'id_number': fid, 
                        'id': fid,
                        'fid_from_user': fid_from_user,
                        'mz': mz, 
                        'rtime': rtime,
                        'pval': p_value,
                        'statistic': statistic,
                        }
                self.ListOfMassFeatures.append( 
                    peak
                    )
            else:
                excluded_list.append( (ii, mz, rtime) )
        
        if excluded_list:
            print( "Excluding %d features out of m/z range %s." %(len(excluded_list), str(MASS_RANGE)) )

        
    def read_from_file(self, inputFile):
        return open(inputFile).read()
    
    def read_from_webform(self, t):
        return t

    def __check_redundant__(self, L):
        redundant = len(L) - len(set(L))
        if redundant > 0:
            print( "Your input file contains %d redundant features." %(redundant) )
        return L

    def read(self):
        '''
        Read input feature lists to ListOfMassFeatures. 
        Row_numbers (rowii+1) are used as primary ID.
        # not using readlines() to avoid problem in processing some Mac files
        '''
        if self.web:
            self.text_to_ListOfMassFeatures(self.paradict['datatext'])
        else:
            self.text_to_ListOfMassFeatures( 
                open(os.path.join(self.paradict['workdir'], self.paradict['infile'])).read() )

        print("Read %d features as reference list." %len(self.ListOfMassFeatures))
    
    
    # more work?
    def determine_significant_list(self, all_feature_list):
        '''
        For single input file format in ver 2. 
        The significant list, input_mzlist, should be a subset of ref_mzlist,
        determined either by user specificed --cutoff,
        or by automated cutoff close to a p-value hotspot,
        in which case, paradict['cutoff'] is updated accordingly.

        '''
        if not self.paradict['cutoff']:
            # automated cutoff
            new = sorted(all_feature_list, key=lambda x: x['pval'])
            
            p_hotspots = [ 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001 ]
            N_hotspots = [ len([x for x in all_feature_list if x['pval'] < pp]) for pp in p_hotspots ]
            
            N_quantile = len(new) / 4
            N_optimum, N_minimum = 300, 30
            chosen = 9999
            for ii in range( len(N_hotspots) ):
                # will get the smallest p as ii increases
                if N_optimum < N_hotspots[ii] < N_quantile:
                    chosen = ii
            
            # if nothing was chosen
            if chosen > 100:
                for ii in range( len(N_hotspots) ):
                    if N_minimum < N_hotspots[ii] < N_quantile:
                        chosen = ii
            
            if chosen > 100:
                N_chosen = int(N_quantile)
                self.paradict['cutoff'] = new[N_chosen+1]['pval']
            else:
                #N_chosen = N_hotspots[chosen]
                
                self.paradict['cutoff'] = p_hotspots[chosen]
        
            print("Automatically choosing (p < %f) as significant cutoff."  %self.paradict['cutoff'])  
        
        # mark MassFeature significant
        for f in self.ListOfMassFeatures:
            if f['pval'] < self.paradict['cutoff']:
                f['is_significant'] = True
        
        self.input_featurelist = [f.row_number for f in self.ListOfMassFeatures if f['is_significant']]
        print("Using %d features (p < %f) as significant list." 
                              %(len(self.input_featurelist), self.paradict['cutoff']))  













class DataMeetModel:
    '''
    changing in v3

    TrioList is no longer relevant

    Changing both metabolic models and expt data to center on neutral mass (formula)

    '''
    def __init__(self, metabolicModel, userData):
        '''
        # from ver 1 to ver 2, major change in .match()
        Trio structure of mapping
        (M.row_number, EmpiricalCompounds, Cpd)
        
        '''
        self.model = metabolicModel
        self.data = userData
        
        # retention time window for grouping, based on fraction of time or ranks
        self.rtime_tolerance = self.data.max_retention_time * RETENTION_TIME_TOLERANCE_FRAC
        self.rtime_tolerance_rank = len(self.data.ListOfMassFeatures) * RETENTION_TIME_TOLERANCE_FRAC
        
        # major data structures
        # web
        if self.data.web:
            wanted_ions = self.data.paradict['wanted_adduct_list']
        # local
        else:
            wanted_ions = wanted_adduct_list[ self.data.paradict['mode'] ]

        self.IonCpdTree = self.__build_cpdindex__(wanted_ions)

        self.rowDict = self.__build_rowindex__( self.data.ListOfMassFeatures )
        self.ListOfEmpiricalCompounds = self.get_ListOfEmpiricalCompounds()
        
        # this is the reference list
        self.mzrows = [M.row_number for M in self.data.ListOfMassFeatures]
        
        self.rowindex_to_EmpiricalCompounds = self.__make_rowindex_to_EmpiricalCompounds__()
        self.Compounds_to_EmpiricalCompounds = self.__index_Compounds_to_EmpiricalCompounds__()
        
        # this is the sig list
        self.significant_features = self.data.input_featurelist
        self.TrioList = self.batch_rowindex_EmpCpd_Cpd( self.significant_features )


    def __build_cpdindex__(self, wanted_ions):
        '''
        indexed Compound list, to speed up m/z matching.
        Limited to MASS_RANGE (default 50 ~ 2000 dalton).
        
        changing from adduct_function to wanted_adduct_list dictionary
        
        wanted_adduct_list['pos_default'] = ['M[1+]', 'M+H[1+]', 'M+2H[2+]', 'M(C13)+H[1+]', 'M(C13)+2H[2+]', 
                    'M+Na[1+]', 'M+H+Na[2+]', 'M+HCOONa[1+]'
                    ],
        
        # 
        >>> metabolicModels['human_model_mfn']['Compounds'].items()[92]
        ('C00217', {'formula': '', 'mw': 147.0532, 'name': 'D-Glutamate; D-Glutamic acid; D-Glutaminic acid; D-2-Aminoglutaric acid',
         'adducts': {'M+2H[2+]': 74.53387646677, 'M+Br81[-]': 227.9695, 'M-H2O+H[1+]': 130.04987646677, 
         'M-C3H4O2+H[1+]': 76.03937646677, 'M-HCOOH+H[1+]': 102.05507646676999, 'M-HCOONa+H[1+]': 80.07307646677, 
         'M+K[1+]': 186.01597646677, 'M+Cl[-]': 182.0221, 'M+Na-2H[-]': 167.02064706646001, 'M-CO2+H[1+]': 104.07067646677, 
         'M+Na[1+]': 170.04247646677, 'M+Br[-]': 225.9715, 'M(S34)-H[-]': 148.04172353323, 'M+H[1+]': 148.06047646677, 
         'M-H4O2+H[1+]': 112.03927646677, 'M(C13)-H[-]': 147.04932353323, 'M(Cl37)-H[-]': 148.04312353323, 'M+HCOONa[1+]': 216.04787646677, 'M(C13)+2H[2+]': 75.03557646677, 'M+HCOOK[1+]': 232.02177646677, 'M-CO+H[1+]': 120.06547646677, 'M+HCOO[-]': 192.050845, 'M(C13)+3H[3+]': 50.359409800103336, 'M(Cl37)+H[1+]': 150.05767646677, 'M-H[-]': 146.04592353323, 'M+ACN-H[-]': 187.07246853323, 'M+Cl37[-]': 184.0191, 'M-H2O-H[-]': 128.03532353322998, 'M(S34)+H[1+]': 150.05627646677002, 'M-HCOOK+H[1+]': 64.09917646677, 'M+3H[3+]': 50.025009800103334, 'M+CH3COO[-]': 206.066495, 'M(C13)+H[1+]': 149.06387646677, 'M[1+]': 147.0532, 'M-NH3+H[1+]': 131.03397646677, 'M+NaCl[1+]': 206.01907646677, 'M+H+Na[2+]': 85.52487646677, 'M+H2O+H[1+]': 166.07107646677002, 'M-H+O[-]': 162.04083353323, 'M+K-2H[-]': 182.99414706646002, 'M-2H[2-]': 72.51932353323001}})
        >>> len(metabolicModels['human_model_mfn']['Compounds'])
        3560



        No more pre-computed adducts




        '''

        IonCpdTree = []
        
        for ii in range(MASS_RANGE[1]+1): 
            IonCpdTree.append([])       #empty lists for anything below MASS_RANGE
            
        # iteritems vs items is contention of efficiency, but there's change btw Python 2 and Python 3...
        for c,d in self.model.Compounds.items():
            if d['mw']:                 #sanity check; bypass mistake in adducts type
                for ion,mass in d['adducts'].items():
                    if ion in wanted_ions and MASS_RANGE[0] < mass < MASS_RANGE[1]:
                        IonCpdTree[ int(mass) ].append( (c, ion, mass) )
                
        # tree: (compoundID, ion, mass), ion=match form; mass is theoretical
        return IonCpdTree


    def __build_rowindex__(self, ListOfMassFeatures):
        '''
        Index list of MassFeatures by row# in input data
        '''
        rowDict = {}
        for M in ListOfMassFeatures: rowDict[M.row_number] = M
        return rowDict


    def __match_all_to_all__(self):
        '''



        change to JMS/khipu based annotation









        Major change of data structure here in version 2.
        In ver 1, matched m/z is stored in each Compound instance.
        Here, we produce mapping dictionaries for
            * mzFeatures to theoretical ions
            * Compounds to mzFeatures
        Then, 
            * EmpiricalCompounds are determined within Compound matched mzFeatures, considering retention time.
        
        
        '''
        self.__match_to_mzFeatures__()
        self.cpd2mzFeatures = self.index_Compounds_to_mzFeatures()
        return self.compound_to_EmpiricalCompounds()
        


    def compound_to_EmpiricalCompounds(self):
        '''
        EmpiricalCompounds are constructed in this function.
        First splitting features matching to same compound by retention time;
        then merging those matched to same m/z features.
        run after self.index_Compounds_to_mzFeatures()
        '''
        ListOfEmpiricalCompounds = []
        for k,v in self.cpd2mzFeatures.items():
            ListOfEmpiricalCompounds += self.__split_Compound__(k, v)      # getting inital instances of EmpiricalCompound
            
        print ("Got %d ListOfEmpiricalCompounds" %len(ListOfEmpiricalCompounds))
        
        # merge compounds that are not distinguished by analytical platform, e.g. isobaric
        return self.__merge_EmpiricalCompounds__( ListOfEmpiricalCompounds )
        


    def __index_Compounds_to_EmpiricalCompounds__(self):
        '''
        Make dict cpd - EmpiricalCompounds
        '''
        mydict = {}
        for E in self.ListOfEmpiricalCompounds:
            for m in E.compounds:
                if m in mydict:
                    mydict[m].append(E)
                else:
                    mydict[m] = [E]
                    
        return mydict
        

    def batch_rowindex_EmpCpd_Cpd(self, list_features):
        '''
        Batch matching from row feature to Ecpds; Use trio data structure, (M.row_number, EmpiricalCompounds, Cpd).
        Will be used to map for both sig list and permutation lists.
        '''
        new = []
        for f in list_features:
            for E in self.rowindex_to_EmpiricalCompounds.get(f, []):
                for cpd in E.compounds:
                    new.append((f, E, cpd))
            
        return new

            
    def get_ListOfEmpiricalCompounds(self):
        '''
        Collect EmpiricalCompounds.
        Initiate EmpCpd attributes here.
        '''
        ListOfEmpiricalCompounds, ii = [], 1
        for EmpCpd in self.__match_all_to_all__():
            EmpCpd.evaluate()
            EmpCpd.EID = 'E' + str(ii)
            EmpCpd.get_mzFeature_of_highest_statistic( self.rowDict )
            ii += 1
            if self.data.paradict['force_primary_ion']:
                if EmpCpd.primary_ion_present:
                    ListOfEmpiricalCompounds.append(EmpCpd)
            else:
                ListOfEmpiricalCompounds.append(EmpCpd)
        
        print ("Got %d final ListOfEmpiricalCompounds" %len(ListOfEmpiricalCompounds))
        return ListOfEmpiricalCompounds

    
    def to_json(self):
        '''
        JSON export to be consumed by downstream functions

        empCpd2Cpds = {empCpd: (), ...,}

        Will update later in accordance to 
        https://github.com/shuzhao-li/metDataModel

        '''

        empCpd2Features, empCpd2Cpds = {}, {}
        for E in self.ListOfEmpiricalCompounds:
            empCpd2Features[E.EID] = E.massfeature_rows
            empCpd2Cpds[E.EID] = E.compounds

        return {
            'metabolic_model': self.model.version,
            'empCpd2Features': empCpd2Features,
            'empCpd2Cpds': empCpd2Cpds,
        }





