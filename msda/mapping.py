import requests
import re
import pandas as pd


def ensp2uid(ensp_id):
    url = "http://grch37.rest.ensembl.org/xrefs/id/"
    query = "%s?content-type=application/json;" % ensp_id
    db = "external_db=Uniprot/SWISSPROT"

    r = requests.get(url + query + db)
    dict = r.json()
    try:
        id = dict[0]['primary_id']
    except IndexError:
        id = 'unknown'
    except KeyError:
        id = 'unknown'
    return id


def name2entrez(name):
    df_map = pd.read_table('resources/gene_name_mapping.txt')
    entrez_id = df_map[df_map['Approved Symbol'] == name][
        'Entrez Gene ID'].values[0]
    return int(entrez_id)


def entrez2name(entrez_id):
    df_map = pd.read_table('resources/gene_name_mapping.txt')
    gene_name = df_map[df_map['Entrez Gene ID'] == entrez_id][
        'Approved Symbol'].values[0]
    return gene_name



def uid2gn(uid):
    uid = uid.split('-')[0]
    url = "http://www.uniprot.org/uniprot/%s.fasta" % uid

    r = requests.get(url)
    try:
        line1 = r.text.split('\n')[0]
        gn = re.search('GN=(.*) PE', line1).group(1)
    except AttributeError:
        gn = 'unknown'
    return gn

def get_uniprot_id(gene_name):
    df_map = pd.read_table('resources/hgnc_mapping.txt')
    try:
        uniprot_id = df_map[df_map['Approved Symbol'] == gene_name][
            'UniProt ID(supplied by UniProt)'].values[0]
    except IndexError:
        uniprot_id = None
    except ValueError:
        uniprot_id = None
    return uniprot_id
