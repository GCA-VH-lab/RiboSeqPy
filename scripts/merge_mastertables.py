#!/usr/bin/env python
#
# Merge master tables of 2 replicates
# Keeps codons covered in both replicates
#
import sys
import argparse
import pandas as pd
import numpy as np
from scipy import stats

parser = argparse.ArgumentParser(description='Merge master tables of 2 replicates')
parser.add_argument('-iR1', type=str, help='Master table for replica 1')
parser.add_argument('-iR2', type=str, help='Master table for replica 2')
parser.add_argument('-out1', type=str, help='Output file name for simply joined master table', default="Master_table_joined_R1R2.txt")
parser.add_argument('-out2', type=str, help='Output file name for updated joined master table', default="Master_table_joined_R1R2_upd.txt")
parser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")
args = parser.parse_args()

if args.verbose:
    print("verbosity turned on")
    
print("\n-iR1  input table replica 1: {}\n-iR2  input table replica 2: {}\n-out1        output_1:          {}\n-out2        output_2:          {}".format(args.iR1, args.iR2, args.out1, args.out2))

usage = "./merge_mastertables.py -iR1 Master_table_MetS2_WTS1_cleaned.txt -iR2 Master_table_MetS4_WTS3_cleaned.txt"

if (args.iR1==None)|(args.iR2==None)|(args.out1==None)|(args.out2==None):
     sys.exit("\n  usage:\n\t{}\n".format(usage))

iR1  = args.iR1
iR1  = args.iR1
out1 = args.out1
out2 = args.out2

# output  merged data file
fo     = args.out1  # just joined
fo_upd = args.out2  # updated - some column names are changed; duplicated columns removed

#R1
df1 = pd.read_csv(args.iR1, sep="\t")
df1.drop('Unnamed: 0', axis=1, inplace=True)
df1['log2(codon_relative_fd)'] = np.log2(df1.codon_relative_fd)
df1['Z-score']= stats.zscore(df1['log2(codon_relative_fd)'])
# R2
df2 = pd.read_csv(args.iR2, sep="\t")
df2.drop('Unnamed: 0', axis=1, inplace=True)
df2['log2(codon_relative_fd)'] = np.log2(df2.codon_relative_fd)
df2['Z-score']= stats.zscore(df2['log2(codon_relative_fd)'])

# create uid for codon gene_id and P-Site codon leftmot position
df1['CidR1']=df1['Gene_id_treated']+"_"+df1['Position_leftmost_1'].astype(str) 
df2['CidR2']=df2['Gene_id_treated']+"_"+df2['Position_leftmost_1'].astype(str)

df1.set_index('CidR1', inplace=True)
df2.set_index('CidR2', inplace=True)

# in the case if index is different between replicas
#index =list(set(df1.index)|set(df2.index))

# sequence+  sequence on a gene if '-' strand then rev comp

df1.loc[df1.Strand=='+','sequence+'] = df1.sequence
df1.loc[df1.Strand=='-','sequence+'] = df1.reverse_complement

df2.loc[df2.Strand=='+','sequence+'] = df2.sequence
df2.loc[df2.Strand=='-','sequence+'] = df2.reverse_complement


col1 = ['Gene_id_treated', 'Chr', 'Strand', 'codon_relative_fd', 'log2(codon_relative_fd)',
        'Z-score', 'sequence+', 'peptide', 'gene_leftmost', 'gene_rightmost', 
        'Position_leftmost_1','u_rpm_1','u_rpm_2', 'codon_E','codon_P','codon_A',
        'codon_relative_rpm_1', 'codon_relative_rpm_2']

df1s = df1[col1]
df2s = df2[col1]

# rename columns
#replica = "R1"
colren_R1 ={i:"{}_R1".format(i) for i in col1}
colren_R2 ={i:"{}_R2".format(i) for i in col1}
# rename columns
df1s.rename(columns=colren_R1, inplace=True)
df2s.rename(columns=colren_R2, inplace=True)
# make new df
data = [df1s, df2s]
df5 = pd.concat(data, axis=1, join='inner') ################### CHEK! common codons!
# remove some duplicated columns
columns_to_remove = ['Gene_id_treated_R2', 'gene_leftmost_R2', 'gene_rightmost_R2',
                     'Position_leftmost_1_R2', 'Chr_R2', 'Strand_R2']

df5.drop(columns_to_remove, axis=1, inplace=True)

#df5.to_csv(fo, sep="\t")
print("DataFrame.shape {}".format(df5.shape))
print("Union of indexes df1 & df2  {}".format(len(set(df1s.index)|set(df2s.index))))
# save outp1
df5.to_csv(fo, sep='\t')
print("Output 1 saved to: {}".format(fo))
#print("Output saved to: {}".format(fo))
print("\nUpdating table ...")

#
# Updating Merged table
#

df = df5.copy()
col2remove = ['codon_E_R2', 'codon_P_R2','codon_A_R2', 'sequence+_R2', 'peptide_R2']

df.drop(col2remove, axis=1, inplace=True)

# list of columns containing "_R1" what are renamed
l = ['Gene_id_treated_R1', 'Chr_R1', 'Strand_R1', 'sequence+_R1',
     'peptide_R1', 'gene_leftmost_R1', 'gene_rightmost_R1','Position_leftmost_1_R1',
     'codon_E_R1', 'codon_P_R1', 'codon_A_R1', 'Position_leftmost_1']

col2rename = {i:i.replace("_R1", "") for i in l}
col2rename['Gene_id_treated_R1'] = 'Gene_id'
col2rename['Position_leftmost_1_R1'] = 'Position_leftmost'
df.rename(columns=col2rename, inplace=True)

colReorder = ['Chr', 'Gene_id',  'Strand', 'codon_relative_fd_R1','codon_relative_fd_R2',
              'u_rpm_1_R2', 'u_rpm_2_R2', 'u_rpm_1_R1', 'u_rpm_2_R1', # remove when not needed
              'log2(codon_relative_fd)_R1','log2(codon_relative_fd)_R2', 'Z-score_R1', 'Z-score_R2', 
              'codon_E','codon_P', 'codon_A', 'sequence+', 'peptide',
              'gene_leftmost', 'gene_rightmost', 'Position_leftmost',
              'codon_relative_rpm_1_R1','codon_relative_rpm_2_R1','codon_relative_rpm_1_R2','codon_relative_rpm_2_R2'
             ]

df = df[colReorder]
print("Updated output 2 saved to: {}".format(fo_upd))
df.to_csv(fo_upd, sep='\t')
print("\nDone!")
