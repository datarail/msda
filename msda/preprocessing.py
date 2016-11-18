import pandas as pd
import re
from mapping import uid2gn

df_map = pd.read_csv('resources/Uniprot_sec_to_prim.csv', sep='\t')

delac_tr = ['C9JYP6']


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
    else:
        print file, "is not a supported filetype"
    columns = [re.sub(" ", "_", str(i)) for i in df.columns]

    if sample_list:
        sample_index = [columns.index(s) for s in columns
                        if "_sn_" in s]
        columns[sample_index[0]: sample_index[-1]+1] = sample_list
    df.columns = columns

    uniprot_id = [re.split('\|', i)[1] for i in df.Protein_Id]
    df.Protein_Id = uniprot_id

    # Find secondary accession numbers in dataframe
    sec_ids = list(set(df.Protein_Id.tolist()).
                   intersection(df_map.Secondary_ID.tolist()))

    # Replace secondary accesion numbers with primary ID's
    for i, pid in enumerate(df.Protein_Id):
        if pid in sec_ids:
            df.Protein_Id.iloc[i] = get_primary_ids(pid)

    # Fix datetime entries in Gene names
    for i, gs in enumerate(df.Gene_Symbol):
        if isinstance(gs, pd.datetime):
            df.Gene_Symbol.ix[i] = uid2gn(df.Protein_Id.ix[i])

    for i, gs in enumerate(df.Gene_Symbol):
        if type(gs) == float:
            df.Gene_Symbol.ix[i] = uid2gn(df.Protein_Id.ix[i])

    # Remove deleted uniprot ids
    df = df[~df.Protein_Id.isin(delac_tr)]

    df['Protein_Id_prefix'] = [pr.strip().split('-')[0]
                               for pr in df.Protein_Id.tolist()]

    return df


def get_primary_ids(secondary_id):
    ind = df_map.Secondary_ID.tolist().index(secondary_id)
    primary_id = df_map.Primary_ID[ind]
    return primary_id
