import pandas as pd
import collections
import re
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from sklearn.cross_validation import Bootstrap
import numpy as np


# http://pir.georgetown.edu/cgi-bin/comp_mw.pl?ids=P53039&seq=&submit=Submit
# http://web.expasy.org/peptide_mass/
# line_ind = [i for i,line in enumerate(lines) if line.startswith('>')]
# ls = [(i,v) for i, v in zip(line_ind[:-1], line_ind[1:])]
# for l in ls:
#     lines2.append(lines[l[0]].strip('\n'))
#     ljoin = ''.join([l.strip('\n') for l in lines[l[0]+1: l[1]]])
#     lines2.append(ljoin)
# with open(filename, 'wb') as f:
#       for line in lines2:
#           f.write('%s\n' % line)
df_mw_chart = pd.read_csv('../data/amino_acid_mw_chart.csv')


def get_id(md_string):
    try:
        id = md_string.strip().split('|')[1]
    except IndexError:
        id = md_string
        print id
    return id


def get_mw(seq):
    counter = collections.Counter(seq)
    mw_seq = 0
    for aa in counter.keys():
        try:
            mw_aa = df_mw_chart['mol_wt(g/mol)'][df_mw_chart['1letter_code']
                                                 == aa].values[0]
            mw_seq += mw_aa * counter[aa]
        except IndexError:
            print aa
    mw_seq -= 18.015 * (len(seq) - 1)
    return mw_seq * 0.001


def get_peptides(seq, amino_acid):
    r_indeces = [m.end()-1 for m in re.finditer(amino_acid, seq)]
    if amino_acid == 'R':
        if r_indeces[-1] < (len(seq)-1):
            r_indeces = [r for r in r_indeces if seq[r+1] != 'P']
        else:
            r_indeces[:-1] = [r for r in r_indeces[:-1] if seq[r+1] != 'P']
    r_indeces.append(len(seq))
    start = 0
    r_peptides = []
    for ind in r_indeces:
        r_peptides.append(seq[start:ind+1])
        start = ind + 1
    return r_peptides


def observable_peptides(seq):
    try:
        r_peptides = get_peptides(seq, 'R')
        rp_peptides = []
        for pep in r_peptides:
            if 'K' in pep:
                k_peptides = pep.split('K')
                k_peptides2 = [p+'K' for p in k_peptides[:-1]]
                rp_peptides += k_peptides2
                rp_peptides.append(k_peptides[-1])
            else:
                rp_peptides.append(pep)
    except IndexError:
        rp_peptides = get_peptides(seq, 'K')
    return rp_peptides


def generate_report(file):
    lines = open(file).readlines()
    uids, length, mw, obs_pep, obs_pep2 = [], [], [], [], []
    for pr in range(0, len(lines), 2):
        if not lines[pr].startswith('>##'):
            pr_id = lines[pr]
            pr_seq = lines[pr+1].strip('\r\n')
            uids.append(get_id(pr_id))
            length.append(len(pr_seq))
            mw.append(get_mw(pr_seq))
            rp_peptides = observable_peptides(pr_seq)
            rp_peptides_subset = [pep for pep in rp_peptides
                                  if 6 < len(pep) < 31]
            obs_pep.append(len(rp_peptides_subset))
            obs_pep2.append(len(list(set(rp_peptides_subset))))
    df = pd.DataFrame(zip(uids, length, mw, obs_pep, obs_pep2),
                      columns=['UniprotID', 'Length',
                               'Molecular_Weight (kDa)',
                               'num_theoretical_peptides', 'num_unique'])
    return df


def compute_ibaq(df, organism='human'):
    ref_file = '../data/%s_proteome_report.csv' % organism
    df_ref = pd.read_csv(ref_file)
    num_theor_peptides, ibaq_list, log10_ibaq = [], [], []
    for protein in df['Protein Id'].tolist():
        uid = protein.strip().split('|')[1]
        num_theor_peptides.append(df_ref[
            df_ref.UniprotID == uid].num_theoretical_peptides.values[0])
    df['num_theoretical_peptides'] = num_theor_peptides
    for id in range(len(df)):
        ibaq = df['default~cq_max_sum'].iloc[id] /\
               df['num_theoretical_peptides'].iloc[id]
        ibaq_list.append(ibaq)
        log10_ibaq.append(np.log10(ibaq))
    df['IBAQ'] = ibaq_list
    df['log10_IBAQ'] = log10_ibaq
    return df


def ups2_regression(ups2_ibaq, ups2_conc):

    regr = LinearRegression()
    regr.fit(ups2_ibaq, ups2_conc)
    # plt.scatter(ups2_ibaq, ups2_conc)
    # m = regr.coef_[0][0]
    # c = regr.intercept_[0]
    # reg_line = [m * x + c for x in ups2_ibaq]
    # plt.plot(ups2_ibaq, reg_line)
    # plt.show()
    return regr


def bootstrap_regression(ups2_ibaq, ups2_conc):
    n_samples = len(ups2_conc)
    N = 10000
    n_dups = N / n_samples
    ups2_conc_sets = []
    ups2_ibaq_sets = []
    for n in range(n_dups):
        ups2_conc_sets += list(ups2_conc)
        ups2_ibaq_sets += list(ups2_ibaq)
    N_set = len(ups2_ibaq_sets)
    bs = Bootstrap(N_set, n_iter=N,
                   train_size=n_samples, random_state=0)
    regr_list = []
    for train_index, _ in bs:
        ibaq = np.array(ups2_ibaq_sets)[train_index]
        conc = np.array(ups2_conc_sets)[train_index]
        regr_list.append(ups2_regression(ibaq, conc))
    return regr_list


def get_mean_sd(regr_list):
    slope = [reg.coef_[0][0] for reg in regr_list]
    slope_mean = np.mean(slope)
    slope_SD = np.std(slope)
    intercept = [reg.intercept_[0] for reg in regr_list]
    intercept_mean = np.mean(intercept)
    intercept_SD = np.std(intercept)
    return slope_mean, slope_SD, intercept_mean, intercept_SD


def calibrate(df, slope, intercept):
    log10_conc = []
    for id in range(len(df)):
        ibaq = df.log10_IBAQ.iloc[id]
        conc = (slope * ibaq) + intercept
        log10_conc.append(conc)
    df['log10_conc'] = log10_conc
    return df
