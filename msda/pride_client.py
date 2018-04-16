import requests

base_url = "http://www.ebi.ac.uk:80/pride/ws/archive/"


def get_peptides(uniprot_id):
    """Main function that returns peptide list for given uniprot id

    Parameters
    ----------
    uniprot_id : str
        The uniprot id to be processed

    Returns
    -------
    peptide_list : list
        list of peptide sequences
    """
    projects = get_project_accession(uniprot_id)
    peptides = []
    for i, project in enumerate(projects):
        print(" processing project %s: %s..." % (i+1, project))
        assay_list = get_assay_accession(uniprot_id, project)
        for assay in assay_list:
            peptides += get_peptide_sequence(assay)
    peptides = list(set(peptides))
    filename = '%s_peptides.txt' % uniprot_id
    with open(filename,  'wb') as f:
        for p in peptides:
            f.write("%s\n" % p)
    return peptides


def get_project_accession(uniprot_id):
    """Send query to PRIDE API for project ids.

    Parameters
    ----------
    uniprot_id : str
        The uniprot id to be processed

    Returns
    -------
    project_list : list
        the list of projects containing peptides for the query uniprot id
    """ 
    params = {'query': uniprot_id,
              'show': '10000',
              'page': '0',
              'order': 'desc',
              'species': '9606'}
    project_list_url = base_url + "project/list"
    r = requests.get(project_list_url, params=params)
    js = r.json()
    project_list = [p['accession'] for p in js['list']
                    if len(p['species']) == 1]
    project_list = list(set(project_list))
    return project_list


def get_assay_accession(uniprot_id, project_accession):
    """Send query to PRIDE API for assay ids

    Parameters
    ----------
    uniprot_id : str
        The uniprot_id to be queried
    project_accession : str
        The project id containing peptides for the above uniprot id

    Returns
    -------
    assay_list : list
         list of two tuple of protein alias and assay id
    """
    query_string = "protein/list/project/%s/protein/%s" % (
        project_accession, uniprot_id)
    assay_list_url = base_url + query_string
    r = requests.get(assay_list_url)
    js = r.json()
    assay_list = [(a['assayAccession'], a['accession']) for a in js['list']]
    return assay_list


def get_peptide_sequence(assay_accession):
    """Send query to PRIDE api for peptide sequence

    Parameters
    ----------
    assay_accesion : tuple
         Two tuple of protein_alias and assay id

    Returns
    -------
    peptide_seq_list : list
         list of peptide sequence
    """
    count_query = "peptide/count/assay/%s" % assay_accession[0]
    count_url = base_url + count_query
    r = requests.get(count_url)
    count = int(r.text)
    pages = count/10000 + 1
    peptide_seq_list = []
    for page in range(pages):
        query_string = "/peptide/list/assay/%s" % assay_accession[0]
        params = {'show': '10000', 'page': page}
        peptide_list_url = base_url + query_string
        r = requests.get(peptide_list_url, params=params)
        js = r.json()
        if js.keys() == ['list']:
            for l in js['list']:
                if assay_accession[1] in l['proteinAccession']:
                    peptide_seq_list.append(l['sequence'])
    return peptide_seq_list
