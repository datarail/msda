import pandas as pd
import os
import numpy as np

path = '../data/test_merge/'
files = os.listdir(path)

df_list = [pd.read_excel("%s%s" % (path, file))
           for file in files]


def remove_zeros(df):
    df = df.rename(columns={'proteinID': 'Protein Id',
                            'geneSymbol': 'gene_symbol'})
    samples = [c for c in df.columns.tolist() if '_sn_sum' in c]
    df2 = df[~(df[samples] == 0).all(axis=1)]
    df3 = df2[~df2['Protein Id'].str.contains('HUMAN_contaminant')]
    return df3


def merge(df_list):
    df_list = [remove_zeros(d) for d in df_list]
    len_df = [len(d) for d in df_list]
    ref_index = len_df.index(max(len_df))
    ref_df = df_list[ref_index]

    ref_samples = [c for c in ref_df.columns.tolist()
                   if '_sn_sum' in c]
    nrm_factor = ref_df[samples].sum().max()
    sample_sum = ref_df[samples].sum().values.tolist()
    nr_list = [nrm_factor / float(s) for s in sample_sum]
    df_nrm_list = [normalize(d, nr_list) for d in df_list]
    df_out = pd.marge([df_nrm_list], ignore_index=True)
    return df_out
    
    
def normalize(df, nrm_list):
    samples =  [c for c in df.columns.tolist() if '_sn_sum' in c]
    df2 = df.copy()
    df2[samples] = df2[samples].multiply(nrm_list)
    return df2
