import pandas as pd
import collections
import re
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
# from sklearn.cross_validation import Bootstrap
import numpy as np
import requests
import os
# import seaborn as sns


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
# df_mw_chart = pd.read_csv('../data/amino_acid_mw_chart.csv')
resource_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'resources')

def get_id(md_string):
    """Uniprot ID extracted from detailed identifiers in fasta file

    Parameters
    -----------
    md_string : str
      input string from fasta file with detailed identifier for each protein

    Returns
    -------
    id : str
      extracted uniprot identifier
    """
    try:
        id = md_string.strip().split('|')[1]
    except IndexError:
        id = md_string
        print(id)
    return id


def get_mw(seq):
    """Calculate molecular weight of input amino acid sequence

    Parameters
    ----------
    seq : str
       sequence of amino acids

    Returns
    -------
    mw_seq : int
      molecular weight of sequence
    """

    counter = collections.Counter(seq)
    mw_seq = 0
    for aa in counter.keys():
        try:
            mw_aa = df_mw_chart['mol_wt(g/mol)'][df_mw_chart['1letter_code']
                                                 == aa].values[0]
            mw_seq += mw_aa * counter[aa]
        except IndexError:
            print(aa)
    mw_seq -= 18.015 * (len(seq) - 1)
    mw_seq = mw_seq * 0.001
    return mw_seq


def get_peptides(seq, amino_acid):
    """in silico digest of sequence given the site at which the enzyme cuts

    Parameters
    ----------
    seq : str
      sequence of amino acids
    amino_acid : str
      one-letter code for site at which enzyme cleaves

    Returns
    -------
    r_peptides: list of strings
       list of petides resulting from in-silico digest
    """
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
    """in silico digest of amino acid sequence by trypsin that cuts at K and R

    Parameters
    ----------
    seq : str
     sequence of amino acids

    Returns
    -------
    rp_peptides : list of strings
       list of tryptic peptides
    """
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
    """proteome fasta file is analyzed to calculate length, molecular weight
    and number of tryptic peptides for all proteins

    Parameters
    ----------
    file : str
       Path to fasta file containing all protein sequences in the organism
      (refer to resources/proteomes)

    Returns
    -------
    df : pandas dataframe
      data table of protein identifier and corresponding mol weight
      and theoretical peptides
    """
    lines = open(file).readlines()
    uids, length, mw, obs_pep, obs_pep2 = [], [], [], [], []
    for pr in range(0, len(lines), 2):
        if (not lines[pr].startswith('>##')) & ('contaminant' not in lines[pr]):
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


def compute_ibaq_1sample(df, organism='human'):
    """IBAQ values computed for total intensities of all proteins

    Parameters
    ----------
    df : pandas dataframe
       proteomics dataset with columns as samples and rows as proteins
    organism : str
       organism (human, rat, or mouse) to calibrate each protein

    Returns
    -------
    df : pandas dataframe
       proteomics dataset normalized by IBAQ
    """
    ref_file = 'resources/%s_proteome_mw_peptides.csv' % organism
    df_ref = pd.read_csv(ref_file)
    num_theor_peptides, ibaq_list, log10_ibaq = [], [], []
    for protein in df['Uniprot_Id'].tolist():
        try:
            uid = protein.strip().split('|')[1]
        except IndexError:
            uid = protein
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


def compute_ibaq_dataset(df, organism='human', samples=None):
    """IBAQ values computed for total intensities of all proteins

    Parameters
    ----------
    df : pandas dataframe
       proteomics dataset with columns as samples and rows as proteins
    organism : str
       organism (human, rat, or mouse) to calibrate each protein

    Returns
    -------
    df : pandas dataframe
       proteomics dataset normalized by IBAQ
    """
    # ref_file = 'resources/%s_proteome_mw_peptides.csv' % organism
    # df_ref = pd.read_csv(ref_file)
    df_ref = pd.read_csv(os.path.join(resource_path,
                                      '%s_proteome_mw_peptides.csv' % organism))
    num_theor_peptides = []
    if samples is None:
        samples = df.columns.tolist()[2:]
    for uid in df['Uniprot_Id'].tolist():
        try:
            num_theor_peptides.append(df_ref[
                df_ref.UniprotID == uid].num_theoretical_peptides.values[0])
        except IndexError:
            r = requests.get('http://www.uniprot.org/uniprot/%s.fasta' % uid)
            seq = r.text.split('\n')[1:]
            seq = ''.join(seq)
            peptides = observable_peptides(seq)
            rp_peptides_subset = [pep for pep in peptides
                                  if 6 < len(pep) < 31]
            num_theor_peptides.append(len(rp_peptides_subset))
    df['num_theoretical_peptides'] = num_theor_peptides
    df2 = df.copy()
    df2 = df2.fillna(0)
    df2[samples] = df2[samples].div(df2['num_theoretical_peptides'],
                                    axis=0)
    
    df2.loc[:, samples] = df2.loc[:, samples].div(df2[samples].sum(axis=0))
    return df2
    # df2[samples] = df2[samples].apply(np.log10)
    # df2[samples] = df2[samples].add(10)
    # df2 = df2.replace([-np.inf], [np.nan])
    # return df2


def ups2_regression(ups2_ibaq, ups2_conc):
    """Linear regression based on concentrations of ups2 standards 
    and their ibaq values

    Parameters
    ----------
    ups2_ibaq : list of floats
       list of log10 ibaq values for UPS2 standards
    ups2_ibaq : list of floats
       list of known concentrations (log10) for UPS2 standards

    Returns
    -------
    regr : tuple
      tuple of slope and intercept values
    """
    ups2_ibaq = np.array(ups2_ibaq).reshape(len(ups2_ibaq), 1)
    ups2_conc = np.array(ups2_conc)
    regr = LinearRegression()
    regr.fit(ups2_ibaq, ups2_conc)
    m = regr.coef_[0]
    c = regr.intercept_
    plt.scatter(ups2_ibaq, ups2_conc)
    reg_line = [m * x + c for x in ups2_ibaq]
    plt.plot(ups2_ibaq, reg_line)
    plt.xlabel('log10 iBAQ')
    plt.ylabel('log10 Conc (fmol)')
    plt.annotate('slope=%.2f, intercept=%.2f' % (m, c),
                 xy=(7, 0))
    plt.show()
    return (m, c)


def bootstrap_regression(ups2_ibaq, ups2_conc):
    """Boostrap regression to get error values on regression parameters

    Parameters
    ----------
    ups2_ibaq : list of floats
       list of ibaq values for UPS2 standards
    ups2_ibaq : list of floats
       list of known concentrations for UPS2 standards

    Returns
    -------
    regr_list : list of tuples
       list of slope and intercept tuples for all run of the boostrap
     """
    ups2_ibaq = np.array(ups2_ibaq).reshape(len(ups2_ibaq), 1)
    ups2_conc = np.array(ups2_conc)
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
    """Summary statistics on bootstrap_regression

    Parameters
    ----------
    regr_list : list of tuples

    Returns
    -------
    tuple : 4-tuple of floats
      mean and SD of slope and intercept
    """
    slope = [reg.coef_[0][0] for reg in regr_list]
    slope_mean = np.mean(slope)
    slope_SD = np.std(slope)
    intercept = [reg.intercept_[0] for reg in regr_list]
    intercept_mean = np.mean(intercept)
    intercept_SD = np.std(intercept)
    return slope_mean, slope_SD, intercept_mean, intercept_SD


def calibrate(df, slope, intercept):
    """Calibrate dataset based on UPS2 standards

    Parameters
    ----------
    df : pandas dataframe
      proteomics dataset normalized by IBAQ
    slope : float
      mean value of slope
    intercept : float
      mean value of intercept

    Returns
    -------
    df : Pandas dataframe
      proteomics data calibrated by ups2 standards to get approximate
      concentrations for all proteins
    """
    log10_conc = []
    for id in range(len(df)):
        ibaq = df.log10_IBAQ.iloc[id]
        conc = (slope * ibaq) + intercept
        log10_conc.append(conc)
    df['log10_conc'] = log10_conc
    return df


def ibaq_comparision(df_ibaq, samples, biomarkers,
                     sample_map=None, plot_name='ibaq_com.png'):

    df1 = df_ibaq[df_ibaq.GENE_NAME.isin(biomarkers)]
    df2 = df1[['GENE_NAME'] + samples]
    df3 = pd.melt(df2, id_vars='GENE_NAME',
                  value_vars=df2.columns.tolist()[1:])
    df3.columns = ['GENE_NAME', 'variable', 'log10(iBAQ)']
    # return df3
    labels = []
    if sample_map is not None:
        for id in range(len(df3)):
            sample = df3['variable'].iloc[id]
            labels.append(sample_map[sample])
        print(len(labels))    
        df3['Label'] = labels        
    #sns.stripplot(x="Gene_Symbol", y="log10(iBAQ)",
    #              data=df3, hue="Label")
    #plt.savefig(plot_name, dpi=600)
    return df3
