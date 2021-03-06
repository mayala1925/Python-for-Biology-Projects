#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 19:02:14 2021

@author: matthewayala
"""


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Reading in Data
gene_expression = pd.read_csv('Files/rpkm_ecdysone_to_ctl.txt', sep = '\t')
response = pd.read_csv('Files/w_ecd_genes.list')
no_response = pd.read_csv('Files/wo_ecd_genes.list')

print(no_response)

# Merging tables to find which genes are controls or not
response_enrich = gene_expression.merge(response, how = 'inner', on = 'Gene_name')
no_response_enrich = gene_expression.merge(no_response, how = 'inner', on = 'Gene_name')

# Adding 'wo': this step was not needed but I thought I would need it.
new_data = []
for i in no_response_enrich['Gene_name']:
    f = 'wo_' + i
    new_data.append(f)

# Creating list in same dimensions of the response and control data frames
yes_ecd = ['Ecd' for y in range(170)]
no_ecd = ['wo_ecd' for y in range(380)]

# This step was not needed
# no_response_enrich['Gene_name'] = new_data

# Adding the "ecd" and 'wo_ecd' columns to each dataframe
no_response_enrich['Ecd'] = no_ecd
response_enrich['Ecd'] = yes_ecd

# Appending response data to the control data
resp_and_no_resp = response_enrich.append(no_response_enrich)
print(resp_and_no_resp)

#Creating the boxplot.
boxxx = sns.boxplot(data = resp_and_no_resp, x = 'Ecd', y = 'Enrichment')
boxxx.set(ylim=(-4, 4))
plt.show()