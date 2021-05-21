#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: matthewayala
"""
import math
import numpy as np
from scipy.stats import fisher_exact

# 1) Writing a function that reads through fasta file and finds all the kmers in the sequences and puts their occurences in a dictionary.

def read_fasta_kmers(file,n):
    fasta_file = open(file,'r')
    sequences = ''
    for i in fasta_file:
        i = i.rstrip()
        if i.startswith('>'):
            i = i.lstrip('>')
            seq_name = i
        else:
            sequence = i
            sequences = sequences + sequence
    
    kmers = {}
    num = len(sequences) - n + 1

    for f in range(num):
        kmer = sequences[f:f+n]
        if kmer not in kmers:
            kmers[kmer] = 0
        kmers[kmer] += 1
    return kmers

enhanced_kmers = read_fasta_kmers('DownstreamIntron.Enhanced.fasta',5)
control_kmers = read_fasta_kmers('DownstreamIntron.Control.fasta',5)

# 2) Finding Log2 Enrichment
# This dictionary turns kmer dictionary with occurences to kmer dictionary with frequencies


def finding_freqs(dictionary):
    total_kmers = 0
    frequencies = []
    for key,value in dictionary.items():
        total_kmers = total_kmers + value
    for k,v in dictionary.items():
        freq = v/total_kmers
        frequencies.append(freq)
    freq_dict = dict(zip(dictionary.keys(),frequencies))
    return freq_dict


enhanced_freqs = finding_freqs(enhanced_kmers)
control_freqs = finding_freqs(control_kmers)

# Turning dictionaries into dataframes -LETTING PANDAS DO ITS THANG

enhanced_df = pd.DataFrame(list(enhanced_freqs.items()),columns = ['Kmers','Freqs'])
control_df = pd.DataFrame(list(control_freqs.items()),columns = ['Kmers','Freqs'])
print(enhanced_df.shape)
print(control_df.shape)
# Merging dataframes together on kmers so enhanced and control frequencies are in the same row.
merged_df = enhanced_df.merge(control_df, on = 'Kmers',suffixes = ['_Enhanced','_Control'])

# Adding Log2_Enrichment column by taking log2 of the enhanced/control freq ratios.
merged_df['Log2_Enrichment'] = np.log2(merged_df['Freqs_Enhanced'] / merged_df['Freqs_Control'])
print(merged_df.head())


# 3) Fisher statistics - This section takes a couple minutesto run 

enhanced_kmer_df = pd.DataFrame(list(enhanced_kmers.items()),columns = ['Kmers','Occurences'])
control_kmer_df = pd.DataFrame(list(control_kmers.items()),columns = ['Kmers','Occurences'])

merged_df_enh = merged_df.merge(enhanced_kmer_df,on = 'Kmers')
merged_df_final = merged_df_enh.merge(control_kmer_df,on = 'Kmers',suffixes = ['_Enhanced','_Control'])

total_enh = merged_df_final['Occurences_Enhanced'].sum()
total_control = merged_df_final['Occurences_Control'].sum()
merged_df_final['Total_Enhanced_kmers'] = total_enh
merged_df_final['Total_Control_kmers'] = total_control

p_values = []
for a,b,c,d in zip(merged_df_final['Occurences_Control'],
                   merged_df_final['Occurences_Enhanced'],
                   merged_df_final['Total_Control_kmers'],
                   merged_df_final['Total_Enhanced_kmers']):
    oddr,pval = fisher_exact([[a,b],[c,d]])
    p_values.append(pval)

# Adding p_values to the dataframe with their respective kmer.

merged_df_final['P_values'] = p_values
print(merged_df_final.head())

# 4) Bonferonni Correction

merged_df_final['Corrected_P_value'] = merged_df_final['P_values'] * 1024

# 5) Sorting kmers by corrected p-values ascending

sorted_df = merged_df_final.sort_values('Corrected_P_value')
print(sorted_df)

