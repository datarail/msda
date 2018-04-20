import pandas as pd
import re
from msda import mapping
import numpy as np
from msda import batch_normalization as bn
import os

resource_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'resources')
df_map = pd.read_csv(os.path.join(resource_path, 'Uniprot_sec_to_prim.csv'),
                     sep='\t')

delac_tr = ['C9JYP6', 'Q7Z469']


def read_dataset(file):
    """Read dataset into a pandas dataframe

    Parameters
    ----------
    file : str

    Returns
    -------
    df : pandas dataframe
    """
    if file.endswith('.xlsx'):
        df = pd.read_excel(file)
    elif file.endswith('.csv'):
        df = pd.read_csv(file)
    elif file.endswith('.tsv'):
        df = pd.read_table(file)
    else:
        raise ValueError("Unknown file extension")
    return df


def verify_column_labels(df, pMS=False):
    columns = df.columns.tolist()
    standard_labels = ['Uniprot_Id', 'Gene_Symbol']
    max_score_cols = [s for s in df.columns.tolist() if 'Max_Score' in s]
    phospho_labels = max_score_cols + ['Site_Position', 'Motif']
    assert set(standard_labels) < set(columns), "Uniport_Id and Gene_Symbol are expected column labels"
    if pMS:
        assert set(phospho_labels) < set(columns), "Max_Score, Site_position and Motif are expected columns in phosphoproteomics datasets"
    return None


def remove_human_contaminants(df):
    df = df[~df.Uniprot_Id.str.contains('HUMAN_contaminant')]
    df.index = range(len(df))
    return df


def remove_reverse_proteins(df):
    df = df[~df.Uniprot_Id.str.contains('##')]
    df.index = range(len(df))
    return df


def get_primary_ids(secondary_id):
    """ replace secondary uniprot ids by primary uniport ids 
    using the lookup table downloaded from Uniprot """
    try:
        primary_id = df_map[df_map.Secondary_ID == secondary_id]['Primary_ID'].values[0]
    except IndexError:
        primary_id = None
    return primary_id


def correct_gene_names(df):
    """ Fix datetime entries in Gene names
    """
    update_symbols = []
    for i, gs in enumerate(df.Gene_Symbol):
        if (not (isinstance(gs, str))) or (':' in gs):
            update_symbols.append(mapping.get_name_from_uniprot(df.Uniprot_Id.iloc[i]))
        else:
            update_symbols.append(gs)
    df.Gene_Symbol = update_symbols
    return df

def correct_uniprot_identifiers(df):
    """ Update secondary uniprot identifier and remove pre- & suffixes
    """
    # remove pre- and suffix from Uniprot Identifiers
    try:
        uniprot_id = [uid.split('|')[1] for uid in df.Uniprot_Id.tolist()]
        df.Uniprot_Id = uniprot_id
    except IndexError:
        pass

    secondary_identifiers = df[df.Uniprot_Id.isin(df_map.Secondary_ID.tolist())][
        'Uniprot_Id'].tolist()
    # Replace secondary accesion numbers with primary ID's
    uids = []
    for pid in df.Uniprot_Id.tolist():
        if pid in secondary_identifiers:
            updated_pid = get_primary_ids(pid)
        else:
            updated_pid = None
        if updated_pid:
            uids.append(updated_pid)
        else:
            uids.append(pid)
    df.Uniprot_Id = uids

    # Remove deleted uniprot ids
    df = df[~df.Uniprot_Id.isin(delac_tr)]
    
    df['Uniprot_Id_prefix'] = [pr.strip().split('-')[0]
                               for pr in df.Uniprot_Id.tolist()]
    return df

    
def preprocess_dataset(file, pMS=False):
    """Dataset is preprocessed to correct for outdated UniProt Identifiers, 
    datetime errors in gene name, and human contaminant proteins.

    Parameters
    ----------
    file : str or pandas dataframe
    pMS : bool
        If True, dataset is a phosphorpoteomics dataset

    Returns
    -------
    df : pandas datframe
    """
    if isinstance(file, pd.DataFrame):
        df = file.copy()
    else:
        df = read_dataset(file)
        
    # Standardize column labels    
    df = df.rename(columns={'Protein_Id': 'Uniprot_Id',
                            'Protein Id': 'Uniprot_Id',
                            'gene_symbol': 'Gene_Symbol',
                            'Uniprot_ID': 'Uniprot_Id'})
    cols = [re.sub(" ", "_", str(i)) for i in df.columns]
    df.columns = cols
    verify_column_labels(df, pMS)
    df = remove_human_contaminants(df)
    df = remove_reverse_proteins(df)
    df = correct_uniprot_identifiers(df)
    df = correct_gene_names(df)
    return df


def noise_filter(df, rep_suffix='rep'):
    diffcut = 1
    filtered_index = []
    for ind in df.index.tolist():
        x1 = df.loc[ind, df.columns.to_series().str.contains(
            '%s1' % rep_suffix).tolist()].values.tolist()
        x2 = df.loc[ind, df.columns.to_series().str.contains(
            '%s2' % rep_suffix).tolist()].values.tolist()
        meandelta = np.mean(np.abs([a-b for a, b in zip(x1, x2)]))
        meanval = [(a+b)/2.0 for a, b in zip(x1, x2)]
        maxdiff = np.max(meanval) - np.min(meanval)
        if (maxdiff - meandelta) > diffcut:
            filtered_index.append(ind)
    df2 = df.loc[filtered_index]
    return df2


def rank_batch_overlap(dfs):
    arr = np.zeros([len(dfs)]*2)
    for i in range(len(dfs)):
        for j in range(len(dfs)):
            uids1 = dfs[i].index.tolist()
            uids2 = dfs[j].index.tolist()
            arr[i, j] = len(set(uids1).intersection(uids2))
    batches = ['batch_%d' % (d+1) for d in range(len(dfs))]
    df = pd.DataFrame(arr, columns=batches, index=batches)
    overlap_score = []
    for b in batches:
        overlap_score.append(df[b].sum() - df[b].loc[b])
    df['overlap'] = overlap_score
    df = df.sort_values('overlap', ascending=False)
    rank = [int(d.split('_')[1]) for d in df.index.tolist()]
    df_rank = [x for y, x in sorted(zip(rank, dfs))]
    return df_rank


def strip_metadata(df, samples):
    """ return dataframe without metadata columns"""
    # df.index = df.Uniprot_Id.tolist()
    df = df[samples]
    df = df.groupby(df.index).sum()
    return df


def rename_labels(df, meta_df, pMS=False):
    """ repalce TMT label by sample name provided in metadata file """
    tmt_labels = meta_df.TMT_label.tolist()
    samples = meta_df.Sample.tolist()
    if pMS:
        df = make_pMS_identifier(df)
        df.index = df['identifier'].tolist()
    else:
        df.index = df.Uniprot_Id.tolist()
    sample_key = {t: s for t, s in zip(tmt_labels, samples)}
    df = df.rename(columns=sample_key)
    df_samples = [s for s in df.columns.tolist() if s in samples]
    if pMS:
        max_score_cols = [s for s in df.columns.tolist() if 'Max_Score' in s]
        metadata_cols = ['Uniprot_Id', 'Gene_Symbol', 'Motif',
                         'Site_Position'] + max_score_cols
        cols = metadata_cols + df_samples
    else:
        cols = ['Gene_Symbol'] + df_samples
    df = df[cols]
    df = df.loc[:, ~df.columns.duplicated()]
    return df, df_samples


def rename_bridge(df, batch_num):
    """ Append batch number to bridge sample label """
    cols = df.columns.tolist()
    cols[-1] = 'Bridge%d' % batch_num
    df.columns = cols
    return df


def merge_batches(dflist, meta_df, pMS=False, norm=False, scale_value=100):
    df_list = []
    for df in dflist:
        # df = pd_import(file)
        if not pMS:
            df = df.drop_duplicates(['Uniprot_Id'])
        df, samples = rename_labels(df, meta_df, pMS)
        df_scaled = scale(df, samples, scale_value)
        # df.index = df.Uniprot_Id.tolist()
        # df = strip_metadata(df, samples)
        df_list.append(df_scaled)

    if norm:
        df_rank = rank_batch_overlap(df_list)
        df_ref = df_rank[0]
        dfA_list = []
        for df in df_list:
            samples = [s for s in df.columns.tolist()
                       if s in meta_df.Sample.tolist()]
            dfA_list.append(bn.normalize_mix(df, df_ref, samples))
        df_rank = rank_batch_overlap(dfA_list)
        dfb_list = []
        for df in dfA_list:
            samples = [s for s in df.columns.tolist()
                       if s in meta_df.Sample.tolist()]
            dfb_list.append(bn.normalize_per_protein(df, df_rank, samples))
        # dfc_list = [rename_bridge(df, num+1) for num, df
        #                 in enumerate(dfb_list)]
        df_merged = pd.concat(dfb_list, axis=1)
    else:
        df_merged = pd.concat(df_list, axis=1)
    df_merged = df_merged.groupby(level=0, axis=1).apply(
          lambda x: x.apply(combine_duplicates, axis=1))
    # df_pms = df_pms[~df_pms.index.str.contains('##')]
    # df_merged = df_merged.convert_objects(convert_numeric=True)
    return df_merged


def make_pMS_identifier(df):
    """ construct unique phosphosite identifier by combining
    gene name, phosphorylated amino acid and position.
    """
    
    identifier = []
    for id in range(len(df)):
        motif = df.Motif.iloc[id].split(';')
        sites = df.Site_Position.iloc[id].split(';')
        uid = df.Uniprot_Id.iloc[id]
        ms = '_'.join(["%s%s" % (m[6], s) for m, s in zip(motif, sites)])
        identifier.append("%s_%s" % (uid, ms))
    df['identifier'] = identifier
    df2 = drop_duplicate_psites(df)
    return df2

def drop_duplicate_psites(df):
    """ remove duplicate phosphosites. Retain site with maximum intensity across samples """
    df2 = df.copy()
    samples = [s for s in df2.columns.tolist() if 'sum' in s]
    df2['max_int'] = df2[samples].sum(axis=1)
    df2 = df2.sort_values('max_int').drop_duplicates(['identifier'],
                                                     keep='last')
    return df2


def combine_duplicates(x):
    return ';'.join(set(x[x.notnull()].astype(str)))


def quantile_normalize(df):
    """ quantile normalize across genes """
    rank_mean = pd.DataFrame(np.sort(df.values, axis=0),
                             index=range(1, len(df)+1),
                             columns=df.columns).mean(axis=1)
    df_quantile = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df_quantile


def get_samples(key_file):
    lines = open(key_file, 'rt').readlines()
    for i, l in enumerate(lines):
        if l.startswith('TMT-126'):
            entry_lines = lines[i:]
    sample_id = {}
    for l in entry_lines:
        id, sample = l.split()
        sample_id[id] = sample
    df = pd.DataFrame(sample_id.items(),
                      columns=['mass_spec_label', 'samples'])
    sample_list = df['samples'].tolist()
    return sample_list


def merge_duplicate_features(df):
    """ Merge duplicate features in dataframe by
    taking their mean values. feature names have to be the index """

    df2 = df.groupby(df.index).mean()
    return df2


def normalize_pMS_by_protein(dfp, dfm, samples, scale_value=100):
    dfp.index = dfp['Uniprot_Id'].tolist()
    dfm.index = dfm['Uniprot_Id'].tolist()
    dfm = dfm.loc[dfp.index.tolist()]
    dfpn = dfp.copy()
    dfpn[samples] = dfpn[samples].div(dfm[samples])
    dfpn = dfpn.replace([np.inf], np.nan).dropna()
    dfpn_scaled = scale(dfpn, samples, scale_value)
    return dfpn_scaled


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
