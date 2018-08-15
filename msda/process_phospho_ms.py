import pandas as pd
import numpy as np
import os


def filter_contaminants_reverse(df):
    """Renames columns and removes proteins that are
    human contaminants or in reverse sequence

    Parameters
    ----------
    df : pandas dataframe
        the mass spec data prior to filtering contaminants

    Returns
    -------
    df4 : pandas dataframe
        mass spec data that has column names properly named
        and contaminant/reverse protein removed.
    """
    df = df.rename(columns={'proteinID': 'Uniprot_Id',
                            'Protein Id': 'Uniprot_Id',
                            'geneSymbol': 'Gene_Symbol',
                            'gene_symbol': 'Gene_Symbol',
                            'Site Position': 'Site_Position',
                            'maxScoreStr': 'max_score',
                            'motifPeptideStr': 'Motif',
                            'Max Score': 'max_score',
                            'sitePosStr': 'Site_Position'})
    samples = [c for c in df.columns.tolist() if '_sn_sum' in c]
    df2 = df[~(df[samples] == 0).all(axis=1)]
    df3 = df2[~df2['Uniprot_Id'].str.contains('HUMAN_contaminant')]
    df4 = df3[~df3['Uniprot_Id'].str.contains('##')]
    df4['Site_Position'] = [str(sp) for sp in df4['Site_Position'].tolist()]
    return df4


def merge(df_list):
    """Normalizes and merge pST and pY single and conmposite datasets

    Parameters
    ----------
    df_list : list of dataframes
       list of pST/ pY  single and composite phospho datasets

    Returns
    -------
    df_out : pandas dataframe
        datastet that has been normalized for variations in total intensity
        values of each protein. pST/pY single/composite datasets are merged.
    """
    # Standardize column names, removeshuman contaminants,
    # and protiens in reverse sequence
    df_list = [filter_contaminants_reverse(d) for d in df_list]

    # Use largest dataset in df_list as refernece dataset
    # to compute normalziation factors
    len_df = [len(d) for d in df_list]
    ref_index = len_df.index(max(len_df))
    ref_df = df_list[ref_index]

    # Divide summed intensities of each sample in reference dataset by
    # The largest summed intensity across all samples to
    # compute normalized intensities
    ref_samples = [c for c in ref_df.columns.tolist()
                   if '_sn_sum' in c]
    nrm_factor = ref_df[ref_samples].sum().max()
    sample_sum = ref_df[ref_samples].sum().values.tolist()
    nr_list = [nrm_factor / float(s) for s in sample_sum]

    # Use normalized intensities across all samples to
    # normalize samples in each dataset
    df_nrm_list = [normalize(d, nr_list) for d in df_list]
    # Concatentate pST/pY single/composite datasets into single dataframe
    df_out = pd.concat(df_nrm_list, ignore_index=True)
    return df_out, nr_list


def normalize(df, nrm_list):
    """ Normalize samples in the phospho dataset by normalized intensities computed
    from the reference (the largest) dataset.

    Parameters
    ----------
    df : pandas dataframe
        pST or pY (single or composite) dataset.
    nrm_list : list of floats
        normalized intensities computed from the reference
        (the largest) dataset.

    Returns
    -------
    df2 : pandas dataframe
        Normalized pST or pY (single or composite) dataset.
    """
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
    dfp.max_score = dfp.max_score.astype(str)
    max_composite = np.max([len(str(s).split(';'))
                            for s in dfp.max_score.tolist()])
    max_score_columns = ['max_score_%d' % (mc+1)
                         for mc in range(max_composite)]
    dfp[max_score_columns] = dfp['max_score'].str.split(';', expand=True)
    dfp[max_score_columns] = dfp[max_score_columns].astype(float)
    dfp2 = dfp[dfp[max_score_columns].gt(max_score_cutoff).any(axis=1)].copy()
    return dfp2


def process(raw_files_path):
    """ Wrapper function that reads in raw excel files from core
    and outputs merged phoshoMS table.

    Parameters
    ----------
    raw_files_path : str
        path to fodler containing excel files from core
    
    Returns
    -------
    dfp : pandas dataframe
      Merged and normalized phosho ms dataframe 
    """
    raw_files = os.listdir(raw_files_path)
    df_list = []
    for rf in raw_files:
        file_path = "%s/%s" % (raw_files_path, rf)
        df_list.append(pd.read_excel(file_path))
    dfp, _  = merge(df_list)
    return dfp
    

def replace_tmt_labels(dfp, dfm):
    """Replace default tmt labels with sample names provided by metadata file
    
    Parameters:
    -----------
    dfp : pandas dataframe
        protein/phosphoprotein dataset

   dfm : pandas dataframe
       metadata table that should contain atleast 2 columns named
       'tmt_label' and 'sample'

    Returns:
    --------
    dfp : pandas dataframe
        protein/phosphoprotein dataset with tmt labels replaced by sample name
    """
    rdict = {t: s for t, s in zip(dfm['tmt_label'].tolist(),
                                  dfm['sample'].tolist())}
    dfp = dfp.rename(columns=rdict).copy()
    return dfp
