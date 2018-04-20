import pandas as pd
import os
import numpy as np

# path = '../data/test_merge/'
# files = os.listdir(path)

# df_list = [pd.read_excel("%s%s" % (path, file))
#           for file in files]


def filter_contaminants_reverse(df):
    df = df.rename(columns={'proteinID': 'Uniprot_Id',
                            'Protein Id': 'Uniprot_Id',
                            'geneSymbol': 'Gene_Symbol',
                            'gene_symbol': 'Gene_Symbol',
                            'Site Position': 'Site_Position',
                            'maxScoreStr': 'Max_Score',
                            'motifPeptideStr': 'Motif',
                            'Max Score': 'Max_Score',
                            'sitePosStr': 'Site_Position'})
    samples = [c for c in df.columns.tolist() if '_sn_sum' in c]
    df2 = df[~(df[samples] == 0).all(axis=1)]
    df3 = df2[~df2['Uniprot_Id'].str.contains('HUMAN_contaminant')]
    df4 = df3[~df3['Uniprot_Id'].str.contains('##')]
    df4['Site_Position'] = [str(sp) for sp in df4['Site_Position'].tolist()]
    return df4


def merge(df_list):
    df_list = [filter_contaminants_reverse(d) for d in df_list]
    len_df = [len(d) for d in df_list]
    ref_index = len_df.index(max(len_df))
    ref_df = df_list[ref_index]

    ref_samples = [c for c in ref_df.columns.tolist()
                   if '_sn_sum' in c]
    nrm_factor = ref_df[ref_samples].sum().max()
    sample_sum = ref_df[ref_samples].sum().values.tolist()
    nr_list = [nrm_factor / float(s) for s in sample_sum]
    df_nrm_list = [normalize(d, nr_list) for d in df_list]
    df_out = pd.concat(df_nrm_list, ignore_index=True)
    return df_out


def normalize(df, nrm_list):
    samples = [c for c in df.columns.tolist() if '_sn_sum' in c]
    df2 = df.copy()
    df2[samples] = df2[samples].multiply(nrm_list)
    return df2


def scale(df, samples, scale_value=100):
    """Return copy of dataframe with a set of columns normalized

    Parameters
    ----------
    df : pandas dataframe
        The data frame to be normalized
    samples : list of str
        Column identifiers to use for indexing the data frame
    scale_value : Optional[float]
        The value to normalize the selected samples to, default: 100

    Returns
    -------
    df2 : pandas dataframe
        The data frame with the given columns normalized to the given value
    """
    df2 = df.copy()
    df2[samples] = df2[samples].div(df2[samples].sum(axis=1),
                                    axis=0)
    df2[samples] = df2[samples].multiply(scale_value)
    return df2


def filter_max_score(dfp, max_score_cutoff=13.0):
    """Filter out peptides that have a max score value less than
    the cutoff value

    Parameters
    ----------
    dfp : pandas dataframe
         the phosphos mass spec data prior to max score filtering
    max_score_cutoff : Optional[float]
         the value below which peptides will be filtered out. Default 13

    Returns
    -------
    dfp2 : pandas dataframe
        phospho mass spec data with peptides that have value above cutoff.
    """
    dfp.Max_Score = dfp.Max_Score.astype(str)
    max_composite = np.max([len(str(s).split(';'))
                            for s in dfp.Max_Score.tolist()])
    max_score_columns = ['max_score_%d' % (mc+1)
                         for mc in range(max_composite)]
    dfp[max_score_columns] = dfp['Max_Score'].str.split(';', expand=True)
    dfp[max_score_columns] = dfp[max_score_columns].astype(float)
    dfp2 = dfp[dfp[max_score_columns].gt(max_score_cutoff).any(axis=1)].copy()
    return dfp2
