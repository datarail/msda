import pandas as pd
import re
from mapping import uid2gn
import numpy as np
import batch_normalization as bn
import os
from msda import process_raw

resource_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'resources')
df_map = pd.read_csv(os.path.join(resource_path, 'Uniprot_sec_to_prim.csv'),
                     sep='\t')

delac_tr = ['C9JYP6', 'Q7Z469']


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


def pd_import(file, sample_list=[]):
    if '.xlsx' in file:
        df = pd.read_excel(file)
    elif '.csv' in file:
        df = pd.read_csv(file)
    elif isinstance(file, pd.DataFrame):
        df = file.copy()
    else:
        print file, "is not a supported filetype"
    df = df.rename(columns={'Protein_Id': 'Uniprot_Id',
                            'Protein Id': 'Uniprot_Id',
                            'gene_symbol': 'Gene_Symbol',
                            'Uniprot_ID': 'Uniprot_Id'})
    cols = [re.sub(" ", "_", str(i)) for i in df.columns]

    if sample_list:
        sample_index = [cols.index(s) for s in cols
                        if "_sn_" in s]
        cols[sample_index[0]: sample_index[-1]+1] = sample_list
    df.columns = cols

    df = df[~df.Uniprot_Id.str.contains('UPSP')]
    try:
        uniprot_id = [re.split('\|', i)[1] for i in df.Uniprot_Id.tolist()]
        df.Uniprot_Id = uniprot_id
    except IndexError:
        pass

    # Find secondary accession numbers in dataframe
    sec_ids = list(set(df.Uniprot_Id.tolist()).
                   intersection(df_map.Secondary_ID.tolist()))

    # Replace secondary accesion numbers with primary ID's
    update_uids = []
    for i, pid in enumerate(df.Uniprot_Id):
        if pid in sec_ids:
            update_uids.append(get_primary_ids(pid))
        else:
            update_uids.append(pid)
    df.Uniprot_id = update_uids

    # Fix datetime entries in Gene names
    update_symbols = []
    for i, gs in enumerate(df.Gene_Symbol):
        if isinstance(gs, pd.datetime):
            update_symbols.append(uid2gn(df.Uniprot_Id.ix[i]))
        elif type(gs) == float:
            update_symbols.append(uid2gn(df.Uniprot_Id.ix[i]))
        elif ':' in gs:
            update_symbols.append(uid2gn(df.Uniprot_Id.ix[i]))
        else:
            update_symbols.append(gs)
    df.Gene_Symbol = update_symbols

    # Remove deleted uniprot ids
    df = df[~df.Uniprot_Id.isin(delac_tr)]

    df['Uniprot_Id_prefix'] = [pr.strip().split('-')[0]
                               for pr in df.Uniprot_Id.tolist()]

    return df


def get_primary_ids(secondary_id):
    """ replace all secondary uniprot ids by primary uniport ids 
    using the lookup table downloaded from Uniprot """
    ind = df_map.Secondary_ID.tolist().index(secondary_id)
    primary_id = df_map.Primary_ID[ind]
    return primary_id


def filter_score(df):
    """ filter phosphosites for which the localization score is below 13 """
    df_split = pd.DataFrame(df.Max_Score.str.split(';').tolist(),
                            index=df.index.tolist()).astype(float)
    df_split = df_split[df_split > 13].dropna(thresh=1)
    df = df.loc[df_split.index.tolist()]
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
        metadata_cols = ['Uniprot_Id', 'Gene_Symbol', 'Motif',
                         'Site_Position', 'Max_Score']
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


def merge_batches(dflist, meta_df, pMS=False, norm=False):
    df_list = []
    for df in dflist:
        # df = pd_import(file)
        if not pMS:
            df = df.drop_duplicates(['Uniprot_Id'])
        df, samples = rename_labels(df, meta_df, pMS)
        df_scaled = process_raw.scale(df, samples)
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
    df_pms = df_pms[~df_pms.index.str.contains('##')]
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



 
