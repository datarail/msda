import pandas as pd
import requests
import numpy as np
import os
import logging

resource_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'resources')
#df_psite = pd.read_table(os.path.join(resource_path,
#                                      'phosphosite_dataset_appended.tsv'))

df_ptm = pd.read_table(os.path.join(resource_path,
                                    'Regulatory_sites_appended.tsv'),
                       error_bad_lines=False)
df_ptm = df_ptm.replace(np.nan, 'nan', regex=True)
#df_ptm = pd.read_csv('resources/ptm_mapping.csv')
df_map = pd.read_csv(os.path.join(resource_path,
                                  'Uniprot_sec_to_prim.csv'), sep='\t')
df_kinase = pd.read_csv(os.path.join(resource_path,
                                     'kinase_substrate_dataset_2016.csv'))


# def get_regulatory_info(uid, msite, mod, organism):
#     """ retrieval of donwnstream regulation of phosohprylation from 
#     phosphosite datasets"""
#     info = {}
#     # df = pd.read_csv('resources/ptm_mapping.csv')
#     df = df_ptm[df_ptm.ORGANISM == organism]
#     site_spec = "%s-%s" % (site, mod)
#     truth_table = (df.ACC_ID == uid) & (df.Site.isin([site_spec]))
#     functions = df.Function[truth_table].values[0]
#     if type(functions) != float:
#         info['Function'] = functions.strip().split(';')
#     else:
#         info['Function'] = 'NA'
#     processes = str(df.Process[truth_table].values[0])
#     if type(processes) != float:
#         info['Process'] = processes.strip().split(';')
#     else:
#         info['Process'] = 'NA'
#     info['PPI'] = df.PPI[truth_table].values[0]
#     info['other'] = df.Other_interactions[truth_table].values[0]
#     pmids = df.PMID[truth_table].values[0]
#     info['evidence'] = pmids.strip().split(';')
#     return info


def get_regulatory_info(motif, organism):
    """ retrieval of donwnstream regulation of phosohprylation from 
    phosphosite datasets"""
    info = {}
    # df = pd.read_csv('resources/ptm_mapping.csv')
    df = df_ptm[df_ptm.ORGANISM == organism]
    df['MOTIF'] = [mtf[1:-1].upper()
                   for mtf in df['SITE_+/-7_AA'].tolist()]
    #id = df.SITE_GRP_ID[df.MOTIF == motif].values[0]
    # info['id'] = get_identifier_by_seq(motif, organism)
    functions = df.ON_FUNCTION[df.MOTIF == motif].values[0]
    info['Function'] = functions.strip().split(';')
    processes = str(df.ON_PROCESS[df.MOTIF == motif].values[0])
    info['Process'] = processes.strip().split(';')
    info['PPI'] = df.ON_PROT_INTERACT[df.MOTIF == motif].values[0]
    info['other'] = df.ON_OTHER_INTERACT[df.MOTIF == motif].values[0]
    pmids = df.PMIDs[df.MOTIF == motif].values[0]
    info['evidence'] = pmids.strip().split(';')
    info['acc_id'] = df.ACC_ID[df.MOTIF == motif].values[0]
    return info


def get_identifier(uid, site, mod, organism):
    df = df_psite[df_psite.ORGANISM == organism]
    # df = pd.read_table('resources/phosphosite_dataset_appended.tsv')
    site_spec = "%s-%s" % (site, mod)
    id = df[(df.MOD_RSD == site_spec) &
            (df.ACC_ID == uid)].SITE_GRP_ID.values[0]
    return id


def get_identifier_by_seq(seq, organism):
    info = {}
    df = df_psite[df_psite.ORGANISM == organism]
    df['MOTIF'] = [mtf[1:-1].upper()
                   for mtf in df['SITE_+/-7_AA'].tolist()]
    try:
        info['id'] = df[df.MOTIF == seq].SITE_GRP_ID.values[0]
        info['uid'] = df[df.MOTIF == seq].ACC_ID.values[0]
    except IndexError:
        info['id'] = 'nan'
        info['uid'] = 'nan'
    return info


def get_kinases(motif, organism='human'):
    df_org = df_kinase[df_kinase.SUB_ORGANISM == organism]
    df_org['MOTIF'] = [mtf[1:-1].upper()
                       for mtf in df_org['SITE_+/-7_AA'].tolist()]
    try:
        kinases = df_org.KINASE[df_org.MOTIF == motif].values[0]
    except IndexError:
        kinases = 'nan'
    return kinases


def get_networkin_kinases(motif):
    df_networkin = pd.read_table('resources/networkin_human_predictions.tsv')
    df_motif = df_networkin[df_networkin.sequnce == motif[1:-1]]
    if not df_motif.empty:
        highest_score = df_motif.networkin_score.max()
        highest_scoring_kinase = df_motif.id[
            df_motif == highest_score].values[0]
    else:
        highest_scoring_kinase == 'nan'
    return highest_scoring_kinase



# def generate_report(input_df, organism):
#     """ report of downstream affects and upstream kinases returned
#     for phosphosite reported in mass spec datasets"""

#     df_input = input_df
#     df_output = df_input.copy()
#     df_output['PSP_idenitifier'] = ['NA']*len(df_output)
#     df_output['Effect_on_function'] = ['NA']*len(df_output)
#     df_output['Effect_on_process'] = ['NA']*len(df_output)
#     df_output['Effect_on_PPI'] = ['NA']*len(df_output)
#     df_output['Evidence'] = ['NA']*len(df_output)

#     for ind in range(len(df_output)):
#         uid = str(df_output.Uniprot_Id.iloc[ind])
#         # aa = str(df_output['amino acid'].iloc[ind])
#         # site_num = df_input['Site Position'].ix[ind]
#         # mod = str(df_input['Mod_Type'].ix[ind])
#         mod = 'p'
#         site = df_output.site.iloc[ind]  # '%s%d' % (aa, site_num)
#         try:
#             id = get_identifier(uid, site, mod, organism)
#             df_output['PSP_idenitifier'].iloc[ind] = id
#             try:
#                 info = get_regulatory_info(uid, site, mod, organism)
#                 df_output['Effect_on_function'].iloc[ind] = info['Function']
#                 df_output['Effect_on_process'].iloc[ind] = info['Process']
#                 df_output['Effect_on_PPI'].iloc[ind] = info['PPI']
#                 df_output['Evidence'].iloc[ind] = info['evidence']
#             except IndexError:
#                 pass
#         except IndexError:
#             pass
#     return df_output

def generate_report(df_input):
    org, psp, func, process, ppi, evidence = [], [], [], [], [], []
    uids = []
    kinases = []
    uid_org = []
    sites = []
    for ind, motif in enumerate(df_input.sequence.tolist()):
        # ids += [df_input.Uniprot_Id.iloc[ind]]*3
        for organism in ['human', 'rat', 'mouse']:
            id_info = get_identifier_by_seq(motif, organism)
            psp.append(id_info['id'])
            uid_org.append(id_info['uid'])
            org.append(organism)
            kinases.append(get_kinases(motif, organism))
            uids.append(df_input.Uniprot_Id.iloc[ind])
            sites.append(df_input.site.iloc[ind])
            try:
                info = get_regulatory_info(motif, organism)               
                func.append(info['Function'])
                process.append(info['Process'])
                ppi.append(info['PPI'])
                evidence.append(info['evidence'])
            except IndexError:
                func.append('nan')
                process.append('nan')
                ppi.append('nan')
                evidence.append('nan')
                
    df_out = pd.DataFrame(zip(uids, sites, org, uid_org, psp, func,
                              process, ppi, evidence, kinases))
    return df_out


def make_input(excel_file):
    """ mass spec datasets are converted to appropriate
    input file (id, site, motif) to query phohsphosite datasets"""
    df_excel = pd.read_excel(excel_file)
    cols = ['Protein Id', 'Site Position', 'sequence']
    df = df_excel[cols].copy()
    df['site'] = ['na']*len(df)
    uids = []
    psites = []
    subseqs = []
    for ind in range(len(df)):
        site_num = df['Site Position'].iloc[ind]
        uniprot_id = df['Protein Id'].iloc[ind].strip().split('|')[1]
        if type(site_num) == unicode:
            sites = site_num.strip().split(';')
        else:
            sites = [site_num]
        uids += [uniprot_id]*len(sites)
        seqs = df['sequence'].iloc[ind].strip().split('#')
        for i in range(len(sites)):
            psites.append("%s%s" % (seqs[i][-1], sites[i]))

    sec_ids = list(set(uids).intersection(df_map.Secondary_ID.tolist()))
    # Replace secondary accesion numbers with primary ID's
    for i, pid in enumerate(uids):
        if pid in sec_ids:
            uids[i] = get_primary_ids(pid)
    for id, site in zip(uids, psites):
        subseqs.append(construct_motif(id, site))
    df_input = pd.DataFrame(zip(uids, psites, subseqs),
                            columns=['Uniprot_Id', 'site', 'sequence'])
    return df_input


def construct_motif(uid, site):

    # logger = logging.getLogger()
    # logger.setLevel(logging.DEBUG)

    # console = logging.StreamHandler()
    # console.setLevel(logging.INFO)
    # formatter = logging.Formatter("%(levelname)s - %(message)s")
    # # tell the handler to use this format
    # console.setFormatter(formatter)
    # logger.addHandler(console)
    
    url = 'http://www.uniprot.org/uniprot/%s.fasta' % uid
    try:
        r = requests.get(url)
        r.raise_for_status()
    except requests.exceptions.HTTPError as err:
        logging.info("%s caused for inquiry protein %s" %  (err, uid))
        raise ValueError('%s returns 404 error' % uid)
    try:
        seq_all = str(r.text)
        site_num = int(site[1:])
        if seq_all == '':
            raise ValueError('Protein %s returns blank string')
        seq = ''.join(seq_all.strip().split('\n')[1:])  
        if len(seq) < site_num:
            logging.info('')
            raise ValueError('Site position %s greater than length of protein %s'
                             % (site_num, uid))
        #if seq_all != '':                        
        st = site_num - 8
        en = site_num + 6
        if st < 0:
            start = 0
            prefix = ''.join(['_']*abs(st))
        else:
            start = st
            prefix = ''
        if en > len(seq):
            end = len(seq) - 1
            suffix = ''.join(['_']*(en - end))
        else:
            end = en
            suffix = ''
        subseq = prefix + seq[start: end+1] + suffix    
    except UnicodeEncodeError:
        print(uid, site)
        subseq = 'obsolote'
    return subseq


def get_primary_ids(secondary_id):
    """ primary ids returned given secondary accesion ID"""
    ind = df_map.Secondary_ID.tolist().index(secondary_id)
    primary_id = df_map.Primary_ID[ind]
    return primary_id


def get_seq(uid):
    """ function that returns seuqnce from uniport give uiprot id"""
    url = 'http://www.uniprot.org/uniprot/%s.fasta' % uid
    r = requests.get(url)
    try:
        seq_all = str(r.text)
        seq = ''.join(seq_all.strip().split('\n')[1:])
    except UnicodeEncodeError:
        print(uid)
        seq = 'unknown'
    return seq


def get_true_site(uid, motif):
    """ function that takes uniprot id and motid and returns
    the phosphosite aminoa cid and position by checking squence 
    uniprot fasta file 
    """
    seq = get_seq(uid)
    if '_' in motif:
        hh = '_'
    elif 'x' in motif:
        hh = 'x'
    aa_position = int(len(motif) / 2.0)
    if motif.startswith(hh):
        sub_motif = motif.strip().split(hh)[-1].upper()
        diff = len(motif) - len(sub_motif)
    elif motif.endswith(hh):
        sub_motif = motif.strip().split(hh)[0].upper()
        diff = 0
    else:
        sub_motif = motif.upper()
        diff = 0
    if seq != 'unknown':
        start_ind = seq.find(sub_motif)
        true_site = start_ind + round(len(motif) / 2.0) - diff
        true_site = "%s%d" % (motif[aa_position - diff], true_site)
    elif seq == 'unknown':
        true_site = 'NA'
    return true_site

#  df_psite_rat[df_psite['SITE_+/-7_AA'].isin(df_inp.sequence.tolist())]


def get_modifications(protein, type='acivity, induced'):
    """ Look up PSP for PTM's that activate/inhibit queried protein """
    df_ptm_human = df_ptm[df_ptm.ORGANISM == 'human']
    df_mods = df_ptm_human[(df_ptm.ACC_ID == protein) &
                           (df_ptm.ON_FUNCTION.str.contains(type))]
    modifications = df_mods.MOD_RSD.tolist()
    return modifications
                                 

