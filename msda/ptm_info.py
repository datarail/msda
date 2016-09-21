import pandas as pd
import requests

df_psite = pd.read_table('resources/phosphosite_dataset_appended.tsv')

df_ptm = pd.read_table('resources/Regulatory_sites_appended.tsv',
                       error_bad_lines=False)
df_map = pd.read_csv('resources/Uniprot_sec_to_prim.csv', sep='\t')
df_kinase = pd.read_csv('resources/kinase_substrate_dataset.csv')


def get_regulatory_info(uid, site, mod, organism):
    """ retrieval of donwnstream regulation of phosohprylation from 
    phosphosite datasets"""
    info = {}
    # df = pd.read_csv('resources/ptm_mapping.csv')
    df = df_ptm[df_ptm.ORGANISM == organism]
    site_spec = "%s-%s" % (site, mod)
    truth_table = (df.Uniprot_ID == uid) & (df.Site.isin([site_spec]))
    functions = df.Function[truth_table].values[0]
    if type(functions) != float:
        info['Function'] = functions.strip().split(';')
    else:
        info['Function'] = 'NA'
    processes = str(df.Process[truth_table].values[0])
    if type(processes) != float:
        info['Process'] = processes.strip().split(';')
    else:
        info['Process'] = 'NA'
    info['PPI'] = df.PPI[truth_table].values[0]
    info['other'] = df.Other_interactions[truth_table].values[0]
    pmids = df.PMID[truth_table].values[0]
    info['evidence'] = pmids.strip().split(';')
    return info


def get_identifier(uid, site, mod, organism):
    df = df_psite[df_psite.ORGANISM == organism]
    # df = pd.read_table('resources/phosphosite_dataset_appended.tsv')
    site_spec = "%s-%s" % (site, mod)
    id = df[(df.MOD_RSD == site_spec) &
            (df.ACC_ID == uid)].SITE_GRP_ID.values[0]
    return id


def get_identifier_by_seq(uid, seq, organism):
    df = df_psite[df_psite.ORGANISM == organism]
    id = df[(df.ACC_ID == uid) &
            (df['SITE_+/-7_AA'])].SITE_GRP_ID.values[0]
    return id


def get_kinase(uid, seq, organism='human'):
    df_org = df_kinase[df_kinase.SUB_ORGANISM == organism]
    try:
        kinases = df_org.KINASE[(df_org.SUB_ACC_ID == uid) &
                                (df_org['SITE_+/-7_AA'] == seq) &
                                (df_org.KIN_ORGANISM == organism)].values[0]
    except IndexError:
        kinases = 'NA'
    return kinases


def generate_report(input_df, organism):
    """ report of downstream affects and upstream kinases returned
    for phosphosite reported in mass spec datasets"""
    
    df_input = input_df
    df_output = df_input.iloc[:, 0:5].copy()
    df_output['PSP_idenitifier'] = ['NA']*len(df_output)
    df_output['Effect_on_function'] = ['NA']*len(df_output)
    df_output['Effect_on_process'] = ['NA']*len(df_output)
    df_output['Effect_on_PPI'] = ['NA']*len(df_output)
    df_output['Evidence'] = ['NA']*len(df_output)

    for ind in range(len(df_output)):
        uid = str(df_output.Uniprot_Id.iloc[ind])
        # aa = str(df_output['amino acid'].iloc[ind])
        # site_num = df_input['Site Position'].ix[ind]
        # mod = str(df_input['Mod_Type'].ix[ind])
        mod = 'p'
        site = df_output.site.iloc[ind]  # '%s%d' % (aa, site_num)
        try:
            id = get_identifier(uid, site, mod, organism)
            df_output['PSP_idenitifier'].iloc[ind] = id
            try:
                info = get_regulatory_info(uid, site, mod, organism)
                df_output['Effect_on_function'].iloc[ind] = info['Function']
                df_output['Effect_on_process'].iloc[ind] = info['Process']
                df_output['Effect_on_PPI'].iloc[ind] = info['PPI']
                df_output['Evidence'].iloc[ind] = info['evidence']
            except IndexError:
                pass
        except IndexError:
            pass
    return df_output


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
    url = 'http://www.uniprot.org/uniprot/%s.fasta' % uid
    r = requests.get(url)
    try:
        seq_all = str(r.text)
        if seq_all != '':        
            seq = ''.join(seq_all.strip().split('\n')[1:])
            site_num = int(site[1:])
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
        else:
            subseq = 'obsolete'
        return subseq
    except UnicodeEncodeError:
        print uid, site
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
        print uid
        seq = 'unknown'
    return seq


def get_true_site(uid, motif):
    """ function that takes uniprot id and motid and returns
    the phosphosite aminoa cid and position by checking squence 
    uniprot fasta file 
    """
    seq = get_seq(uid)
    if motif.startswith('_'):
        motif = motif.strip().split('_')[-1].upper()
        diff = 15 - len(motif)
    elif motif.endswith('_'):
        motif = motif.strip().split('_')[0].upper()
        diff = 0
    else:
        motif = motif.upper()
        diff = 0
    if seq != 'unknown':
        start_ind = seq.find(motif)
        true_site = start_ind + 8 - diff
        true_site = "%s%d" % (motif[7 - diff], true_site)
    elif seq == 'unknown':
        true_site = 'NA'
    return true_site

#  df_psite_rat[df_psite['SITE_+/-7_AA'].isin(df_inp.sequence.tolist())]

