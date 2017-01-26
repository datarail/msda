import requests
import re


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
