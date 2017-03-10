import requests
import pandas as pd
import subprocess


def get_pathsbetween(source_list, filter=False):
    base_url = 'http://www.pathwaycommons.org/pc2/graph'
    source_str = ','.join(source_list)
    params = {'source': source_str,
              'kind': 'PATHSBETWEEN',
              'format': 'EXTENDED_BINARY_SIF'}
    r = requests.get(base_url, params=params)
    with open('temp_pc_interactions.txt', 'wb') as f:
        f.write(r.text)
    df = pd.read_table('temp_pc_interactions.txt')
    if filter:
        df = df[(df.PARTICIPANT_A.isin(source_list)) & (
            df.PARTICIPANT_B.isin(source_list))]
    return df


def make_network_plot(weights_file, network_file, figure_file, subsets=None):
    df = pd.read_csv(weights_file, sep='\t')
    source_list = df['0'].tolist()
    df = get_pathsbetween(source_list, filter=True)
    df.to_csv(network_file, index=False)
    df2 = df[df.columns[:3]]
    if subsets:
        df2 = df2[df2.INTERACTION_TYPE.isin(subsets)]
    df2.to_csv('temp_pc_file.txt', sep='\t', index=False)
    subprocess.call(['Rscript', 'netcontext.R', '-n', 'temp_pc_file.txt',
                     '-w',  weights_file, '-o', network_file])
    return df2

