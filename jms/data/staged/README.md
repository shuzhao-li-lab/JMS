# Staged intermediary files

## json files GEM compound table
- Notice that these json files don't have neutral mass at this stage
- Bigg even don't have charged formula.
```
├── modelSEED.json
├── Bigg.json
└── vmh.json
```

## serum metabolites

- AS generated serum_metabolites_convalues_unique_editsAS.csv using notebooks @ https://github.com/shuzhao-li/JMS/tree/main/notebooks/RetreiveSerumMetabolitesConcValues 
  - which needs to be converted to dictionary, with meta data on units. 

```
├── serum_metabolites_convalues_unique_editsAS.csv
```

## signature mass json
```
├── ModelSeed_universal_signatures_dict.json
```

## Qstd annotation table
- Qstd annotation table that has verified by searching MS2 of Qstd samples on collected authentic standard MS2 spectra
```
├── Qstd_annot_HILICpos_2022-04-18_MG.json
├── Qstd_annot_RPneg_2022-04-18_MG.json
```

## Authentic standard library
- Authentic standard library tabular tables were wrapped into Empirical compound json format.
  - `combn_1strnd` tables contain both MG and SJ's authentic standards. 
```
├── Auth_Std_library_HILICneg_2022-04-18_Shujian\ Zheng|Minghao\ Gong.json
├── Auth_Std_library_Rppos_2022-04-18_Shujian\ Zheng|Minghao\ Gong.json
├── combn_1strnd_HILICpos_2022-04-19_Shujian\ Zheng|Minghao\ Gong.json
├── combn_1strnd_Rpneg_2022-04-19_Shujian\ Zheng|Minghao\ Gong.json
```