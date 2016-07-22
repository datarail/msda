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
