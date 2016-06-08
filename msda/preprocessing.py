import pandas as pd
import re
import numpy as np
import requests
import StringIO


mapping_url = 'http://www.uniprot.org/mapping/'

df_map = pd.read_csv('resources/Uniprot_sec_to_prim.csv', sep='\t')


def pd_import(excel_file, sample_list):
    df = pd.read_excel(excel_file)
    columns = [re.sub(" ", "_", str(i)) for i in df.columns]
    
    sample_index = [columns.index(s) for s in columns
                    if "_sn_scaled" in s]
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
            df.Protein_Id[i] = get_primary_ids(pid)

    # Fix datetime entries in Gene names
    for i, gs in enumerate(df.Gene_Symbol):
        if isinstance(gs, pd.datetime):
            df.Gene_Symbol.ix[i] = get_gene_names(df.Protein_Id.ix[i])
    
    return df


def get_gene_names(uid):

    uid = re.split('\-', uid)[0]

    # query_string = " ".join(uid_list)

    params = {
        'from': 'ACC',
        'to': 'GENENAME',
        'format': 'tab',
        'query': uid
        }

    r = requests.get(mapping_url, params)
    r.raise_for_status()
    mp = r.text

    cmp = StringIO.StringIO(mp)

    df = pd.read_csv(cmp, sep='\t')
    name = re.split('\_', df.To[0])[0]
    return name


def get_primary_ids(secondary_id):
    ind = df_map.Secondary_ID.tolist().index(secondary_id)
    primary_id = df_map.Primary_ID[ind]
    return primary_id



    
    
