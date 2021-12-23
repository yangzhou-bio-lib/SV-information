#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import re


AF_CSE_Fst = pd.read_table('./AF_CSE_Fst.weir.fst')
di_data = AF_CSE_Fst
di_data = di_data.drop('WEIR_AND_COCKERHAM_FST',axis=1)
di_data


fst_list = ["AF_CSE_Fst.weir.fst","AF_HW_Fst.weir.fst","AF_INA_Fst.weir.fst",
            "AF_SC_Fst.weir.fst","AF_WE_Fst.weir.fst","CSE_HW_Fst.weir.fst",
            "CSE_INA_Fst.weir.fst","CSE_SC_Fst.weir.fst","CSE_WE_Fst.weir.fst","HW_INA_Fst.weir.fst","HW_SC_Fst.weir.fst",
            "HW_WE_Fst.weir.fst","INA_SC_Fst.weir.fst","INA_WE_Fst.weir.fst","SC_WE_Fst.weir.fst"]

for fst in fst_list:
    name_two = re.sub("_Fst\.weir\.fst","",fst)
    fst_data = pd.read_table(fst,names=["CHROM","POS",name_two])
    fst_data = fst_data[1:]
    di_data.insert(2,name_two,fst_data[name_two].values)

di_data_replace=di_data.replace(np.nan,0)
di_data_replace_revise = di_data_replace.astype(str).replace(r"-.*",0,regex=True)
di_data_replace_revise = di_data_replace_revise.astype(float)
di_result = AF_CSE_Fst.drop('WEIR_AND_COCKERHAM_FST',axis=1)

for di in fst_list:
#     计算di，并入表格
    name_two = re.sub("_Fst\.weir\.fst","",di)
    fst_data = ((di_data_replace_revise[name_two] - di_data_replace_revise[name_two].mean())/di_data_replace_revise[name_two].std())
    di_result.insert(2,name_two,fst_data.values)

di_result


di = AF_CSE_Fst.drop('WEIR_AND_COCKERHAM_FST',axis=1)

fst_list_chg = ["AF_CSE","AF_HW","AF_INA","AF_SC","AF_WE","CSE_HW","CSE_INA","CSE_SC","CSE_WE","HW_INA","HW_SC","HW_WE","INA_SC","INA_WE","SC_WE"]

for iterm in ['AF', 'CSE', 'HW', 'INA', 'SC', 'WE']:
    iterm_di = pd.Series([0]*14942)
    for two_population in fst_list_chg:
        if iterm in two_population:
            iterm_di = iterm_di + di_result[two_population]
    di.insert(2, iterm, iterm_di)
di

di.to_csv('./sample_339_14942_di_6_population.csv')
