import requests
import re
import pandas as pd
import os
import numpy as np

resource_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'resources')
df_map = pd.read_table(os.path.join(resource_path, 'hgnc_mapping.txt'))


def get_uniprot_from_ensembl(ensp_id):
    """Return Uniprot ID given Ensembl Protein ID

    Parameters
    ----------
    ensp_id : str
        Ensembl Protein ID

    Returns
    -------
    uniprot_id : str
       Uniprot_Id
    """
    url = "http://grch37.rest.ensembl.org/xrefs/id/"
    query = "%s?content-type=application/json;" % ensp_id
    db = "external_db=Uniprot/SWISSPROT"

    r = requests.get(url + query + db)
    dict = r.json()
    try:
        uniprot_id = dict[0]['primary_id']
    except IndexError:
        uniprot_id = None
    except KeyError:
        uniprot_id = None
    return uniprot_id


def get_entrez_from_name(name):
    """Return Entrez Id given Gene Name

    Parameters
    ----------
    name : str
        Gene Name/Symbol

    Returns
    -------
    id : int
       Entrez ID
    """
    try:
        entrez_id = int(df_map[df_map['Approved Symbol'] == name][
            'Entrez Gene ID'].values[0])
    except IndexError:
        entrez_id = None
    return entrez_id


def get_name_from_entrez(entrez_id):
    """Return Gene name given Entrez ID

    Parameters
    ----------
    entrez_id : int
        Entrez ID

    Returns
    -------
    id : str
       Gene Name/Symbol
    """
    try:
        gene_name = df_map[df_map['Entrez Gene ID'] == entrez_id][
            'Approved Symbol'].values[0]
    except IndexError:
        gene_name = None
    return gene_name



def get_name_from_uniprot(uniprot_id):
    """Return Gene name given Uniport Id

    Parameters
    ----------
    uid : str
       Uniprot ID

    Returns
    -------
    gene_name : str
       Gene Name/ Symbol
    """
    uniprot_id= uniprot_id.split('-')[0]
    url = "http://www.uniprot.org/uniprot/%s.fasta" % uniprot_id

    r = requests.get(url)
    try:
        line1 = r.text.split('\n')[0]
        gene_name = re.search('GN=(.*) PE', line1).group(1)
    except AttributeError:
        gene_name = None
    return gene_name


def get_uniprot_from_name(gene_name):
    """Return Uniprot Id given Gene Name

    Parameters
    ----------
    gene_name : str
       Gene Name/Symbol

    Returns
    -------
    uniprot_id : str
       Uniprot ID
    """ 
    try:
        uniprot_id = df_map[df_map['Approved Symbol'] == gene_name][
            'UniProt ID(supplied by UniProt)'].values[0]
    except IndexError:
        uniprot_id = None
    except ValueError:
        uniprot_id = None
    return uniprot_id


def get_name_from_ensembl(ensembl_id):
    """Return Gene name given Entrez ID

    Parameters
    ----------
    ensembl_id : str
        Ensembl ID

    Returns
    -------
    id : str
       Gene Name/Symbol
    """
    try:
        gene_name = df_map[df_map['Ensembl Gene ID'] == ensembl_id][
            'Approved Symbol'].values[0]
    except IndexError:
        gene_name = None
    return gene_name


def get_name_from_synonyms(synonym):
    """Return Gene name given Synonym

    Parameters
    ----------
    synonym : str
        Synonym

    Returns
    -------
    id : str
       Gene Name/Symbol
    """
    df_map2 = df_map.replace([np.nan], 'nan')
    try:
        gene_name = df_map2[df_map2['Synonyms'].str.contains(
            synonym)]['Approved Symbol'].values[0]
    except IndexError:
        gene_name = None
    return gene_name


def get_entrez_from_ensembl(ensembl_id):
    """ Return Gene name given Entrez ID 
    
    Parameter
    ---------
    ensembl_id: str
        Ensembl ID
    
    Return
    ------
    entrez_id: string
       Entrez ID
    """
    try:
        entrez_id = df_map[df_map['Ensembl Gene ID'] == ensembl_id][
             'Entrez Gene ID'].values[0]
        entrez_id = int(entrez_id)
    except IndexError:
        entrez_id = None
    except ValueError:
        entrez_id = None
    return entrez_id
