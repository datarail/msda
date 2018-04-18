import json
import requests
import pandas as pd


def post_gene_set(gene_set, description=''):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    genes_str = '\n'.join(gene_set)
    description = description
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    data['job_name'] = description
    return data


def get_enrichment_result(data, library='KEGG_2015'):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
    query_string = '?userListId=%s&filename=%s&backgroundType=%s'
    user_list_id = data['userListId']
    filename = data['job_name']
    gene_set_library = library

    url = ENRICHR_URL + query_string % (user_list_id,
                                        filename,
                                        gene_set_library)
    response = requests.get(url, stream=True)

    response_lines = response.text.split('\n')
    header = response_lines[0].split('\t')
    body = [x.split('\t') for x in response_lines[1:]]
    try:
        df = pd.DataFrame(body, columns=header)
        return df
    except AssertionError:
        print(data['job_name'])
     #return df
