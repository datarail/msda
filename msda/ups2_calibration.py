import pandas as pd
import collections
import re

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
    id = md_string.strip().split('|')[1]
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
        r_indeces = [r for r in r_indeces if seq[r+1] != 'P']
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
        rp_peptides = get_peptides(seq, 'P')
    return rp_peptides


def generate_report(file):
    lines = open(file).readlines()
    uids, length, mw, obs_pep, obs_pep2 = [], [], [], [], []
    for pr in range(0, len(lines), 2):
        pr_id = lines[pr]
        pr_seq = lines[pr+1].strip('\r\n')
        uids.append(get_id(pr_id))
        length.append(len(pr_seq))
        mw.append(get_mw(pr_seq))
        rp_peptides = observable_peptides(pr_seq)
        rp_peptides_subset = [pep for pep in rp_peptides if 6 < len(pep) < 31]
        obs_pep.append(len(rp_peptides_subset))
        obs_pep2.append(len(list(set(rp_peptides_subset))))
    df = pd.DataFrame(zip(uids, length, mw, obs_pep, obs_pep2),
                      columns=['UniprotID', 'Length',
                               'Molecular_Weight (kDa)',
                               'num_theoretical_peptides', 'num_unique'])
    return df
