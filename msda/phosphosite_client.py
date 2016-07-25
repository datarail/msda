import requests
import pandas as pd

base_url = 'http://www.phosphosite.org/bulkSequenceSearchSubmitAction.action'


def get_ptms(query_string):
    r = requests.post(base_url, files={'sequenceFile': ('tmp_file.txt',
                                                        query_string,
                                                        'text/plain')}).text

    r_list = r.strip().split('\n')

    sequences = [line.strip().split('\t')[1] for line in r_list]
    proteins = [line.strip().split('\t')[2] for line in r_list]
    sites = [line.strip().split('\t')[5] for line in r_list]
    sites = [x if x != '' else None for x in sites]

    df_ptm = pd.DataFrame(zip(sequences, proteins, sites),
                          columns=[sequences[0], proteins[0], sites[0]])

    return df_ptm
