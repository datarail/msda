import re
import pandas as pd


def prune_list(peptide_list):
    tryptic_petides = [p for p in peptide_list if _is_tryptic(p)]
    l1, l2, l3, l4, l5 = [], [], [], [], []
    for pep in tryptic_petides:
        score = score_(pep)
        l1.append(score[0])
        l2.append(score[1])
        l3.append(score[2])
        l4.append(score[3])
        l5.append(score[4])
    df = pd.DataFrame(zip(tryptic_petides, l1, l2, l3, l4, l5),
                      columns=['sequence', 'starts_with_EorDorQ',
                               'starts_with_KKorRRorRK',
                               'ends_with_KEorKDorRDorRE',
                               'length_crit', 'score'])
    return df


def _is_tryptic(peptide_seq):
    tr = not verify_cm(peptide_seq) and verify_kr_end(
        peptide_seq) and not verify_kr_inner(peptide_seq)
    return tr


def verify_cm(peptide_seq):
    cm = bool(re.search('C', peptide_seq)) or bool(
        re.search('M', peptide_seq))
    return cm


def verify_kr_end(peptide_seq):
    kr_end = bool(re.search('K$', peptide_seq)) or bool(
        re.search('R$', peptide_seq))
    return kr_end


def verify_kr_inner(peptide_seq):
    kr_inner = bool(re.search('K', peptide_seq[:-1])) or bool(
        re.search('R', peptide_seq[:-1]))
    return kr_inner


def score_(peptide_seq):
    edq, kr, ked, length_crit, score = False, False, False, False,  0
    if bool(re.search("^[E,D,Q]", peptide_seq)):
        edq = True
        score -= 1
    if bool(re.search("^KK", peptide_seq)) or bool(
            re.search("^RR", peptide_seq)) or bool(
                re.search("^RK", peptide_seq)):
        kr = True
        score -= 1
    if bool(re.search("KE$", peptide_seq)) or bool(
            re.search("KD$", peptide_seq)) or bool(
                re.search("RD$", peptide_seq)) or bool(
                    re.search("RE$", peptide_seq)):
        ked = True
        score -= 1
    if not 5 <= len(peptide_seq) < 35:
        length_crit = True
        score -= 1
    return edq, kr, ked, length_crit, score


def check_ptm_redundancy(df, filename):
    df_ptm = pd.read_table(filename, index_col=False)
    df_ptm['query_sequence'] = [seq.upper()
                                for seq in df_ptm.Sequence.tolist()]
    site_list = []
    redundancy_list = []

    for seq in df.sequence.tolist():
        red = df_ptm.Protein[df_ptm.query_sequence == seq].tolist()
        redundancy_list.append([pr for pr in red if type(pr) == str])
        sites = df_ptm.Site[df_ptm.query_sequence == seq].tolist()
        site_list.append([s for s in sites if type(s) == str])

    for i, l in enumerate(redundancy_list):
        if not l:
            redundancy_list[i] = 'unknown'
            site_list[i] = 'unknown'

    for i, ptm in enumerate(site_list):
        if not ptm:
            site_list[i] = 'no PTM reported'

    df_update = df.dopy()
    df_update['redundancy'] = redundancy_list
    df_update['PTM_sites'] = site_list
        
    return df_update
        

