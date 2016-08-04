import pandas as pd
import requests
import StringIO
import xml.etree.ElementTree as ET
import re


def get_kinase_set(df):

    # From phosphosite plus
    # =====================
    # df_ks = pd.read_csv('resources/Kinase_Substrate_Dataset.csv',
    #                     sep='\t', header=2)
    # kinase_list = list(set(df_ks.KIN_ACC_ID.tolist()))

    # From Kinome website 2007 paper
    # ==============================
    # df_kin = pd.read_excel('resources/Kincat_Hsap.xls')
    # df_s = df_kin[df_kin['Pseudogene?'] == 'N']
    # entrez_list = df_s['Entrez_GENEID']
    # entrez_list = [str(int(i)) for i in entrez_list if not np.isnan(i)]
    # uid_list = [get_uniprot_id(id) for id in entrez_list]
    # uid_list = [i for i in uid_list if i not in ['not found']]
    # joint_list = list(set(kinase_list) | set(uid_list))

    df_ks = pd.read_excel('resources/kinome_tree.xls')
    kinase_list = list(set(df_ks.UniprotID.tolist()))
    df_subset = df.loc[df.Protein_Id.isin(kinase_list)]
    return df_subset


def get_tf_set(df):
    tr_file = open('resources/transcription_factors_uniprot.txt').readlines()
    tr_list = [i.split('\n')[0] for i in tr_file]
    df_subset = df.loc[df.Protein_Id.isin(tr_list)]
    return df_subset


def get_pfam_domain_acc(uniprot_id):
    url = 'http://pfam.xfam.org/protein/%s?output=xml' % uniprot_id
    r = requests.get(url)
    r.raise_for_status()
    xml = r.text
    root = ET.fromstring(xml)

    for child in root[0]:
        if 'matches' in child.tag:
            pfam_matches = child

    if pfam_matches:
        domain_acc = []
        for match in pfam_matches:
            domain_acc.append(match.attrib['accession'])
    return domain_acc


def get_subset(list, domain_list):

    # protein_list = df.Protein_Id.tolist()
    protein_subset = []

    for protein in list:
        protein = re.split('\-', protein)[0]
        domain_acc = get_pfam_domain_acc(protein)
        common_list = set(domain_list).intersection(domain_acc)
        if common_list:
            protein_subset.append(protein)
    return protein_subset


def get_apoptotic_subset(df):
    df_tomahaq = pd.read_excel('resources/apoptotic_tomahaq_targets.xlsx')
    uids = df_tomahaq['Uniprot ID'].tolist()
    uids.remove('Q96PG8, Q9BXH1')
    uids.append('Q96PG8')
    uids.append('Q9BHX1')
    df_subset = df.loc[df.Protein_Id_tr.isin(uids)]
    return df_subset


def get_membrane_subset(df):
    df_membrane = pd.read_table('resources/processed__BrCa_Sets1_12_AtoD_EtoH_AtoD2_EtoH2_DrugsA_B_Tryps_140806_Quant.tsv')
    uids = df_membrane['Uniprot-l'].tolist()
    df_subset = df.loc[df.Protein_Id.isin(uids)]
    return df_subset


def get_uniprot_id(id, db):
    ''' mapping for proteins in kinome.org dataset '''

    mapping_url = 'http://www.uniprot.org/mapping/'

    params = {
        'from': db,
        'to': 'ACC',
        'format': 'tab',
        'query': id
        }

    r = requests.get(mapping_url, params)
    r.raise_for_status()
    mp = r.text

    cmp = StringIO.StringIO(mp)
    df = pd.read_csv(cmp, sep='\t')
    try:
        name = re.split('\_', df.To[0])[0]
    except IndexError:
        # print refseq
        name = 'not found'
    return name
