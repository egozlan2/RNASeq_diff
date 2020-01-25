#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import scipy.stats as stats
import numpy as np


#path to RNAseq file (make sure if patient_ids are duplicated, expression values are averaged & only one remains per barcode):
rna = "/PATH_TO_RNASEQ_FILE"
data = pd.read_csv(rna, index_col = "Hugo_Symbol")

#select your file with groups of interest here:
groups = "/PATH_TO_GROUPS"
data2 = pd.read_excel(groups, sheet_name = "set1")
data3 = pd.read_excel(groups, sheet_name = "set2")

#create list of barcodes to compare:
set1_group1 = data2.iloc[1:,0].dropna().values.tolist()
set1_group2 = data2.iloc[1:,1].dropna().values.tolist()
set2_group1 = data3.iloc[1:,0].dropna().values.tolist()
set2_group2 = data3.iloc[1:,1].dropna().values.tolist()

s1g1_rna = data.filter(items = set1_group1)
s1g2_rna = data.filter(items = set1_group2)

s2g1_rna = data.filter(items = set2_group1)
s2g2_rna = data.filter(items = set2_group2)


# statistics
s1g1_mean = s1g1_rna.mean(axis =1, skipna = True)
s1g2_mean = s1g2_rna.mean(axis =1, skipna = True)

s2g1_mean = s2g1_rna.mean(axis =1, skipna = True)
s2g2_mean = s2g2_rna.mean(axis =1, skipna = True)

s1_log_ratio = np.log10(s1g1_mean/s1g2_mean)
s2_log_ratio = np.log10(s2g1_mean/s2g2_mean)


#t_test
s1_ttest = stats.ttest_ind(s1g1_rna,s1g2_rna, axis = 1)
s1_df_stats = pd.DataFrame(s1_ttest).transpose().drop([0], axis =1)

s2_ttest = stats.ttest_ind(s2g1_rna,s2g2_rna, axis = 1)
s2_df_stats = pd.DataFrame(s2_ttest).transpose().drop([0], axis =1)


with pd.ExcelWriter('/Users/Goozman/Desktop/SIDE_PROJECT/Mary/Mary_RNA_test.xlsx') as writer: 

    s1g1_mean.to_excel(writer, sheet_name = "stats_set1", startcol= 0,index = True, header = ["Mean Set 1, Group 1"])
    s1g2_mean.to_excel(writer, sheet_name = "stats_set1", startcol= 2,index = False,header = ["Mean Set 1, Group 2"])    
    s1_df_stats.to_excel(writer, sheet_name = "stats_set1", startcol= 3, index = False, header = ["P-value"])
    s1_log_ratio.to_excel(writer, sheet_name = "stats_set1", startcol= 4,index = False,header = ["Log Ratio"])
    
    s2g1_mean.to_excel(writer, sheet_name = "stats_set2", startcol= 0,index = True, header = ["Mean Set 2, Group 1"])
    s2g2_mean.to_excel(writer, sheet_name = "stats_set2", startcol= 2,index = False,header = ["Mean Set 2, Group 2"])    
    s2_df_stats.to_excel(writer, sheet_name = "stats_set2", startcol= 3, index = False, header = ["P-value"])
    s2_log_ratio.to_excel(writer, sheet_name = "stats_set2", startcol= 4,index = False,header = ["Log Ratio"])
