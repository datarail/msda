import requests
import io
import pandas as pd
import re
import numpy as np


def get_mutation_data(gene_list, cancer_subtypes):
    """retrieve case-level data on mutations for given list 
    of genes and cell lines
    """
    
    base_url = 'http://www.cbioportal.org/webservice.do'

    genes = ' '.join(gene_list)
    subtypes = ' '.join(['%s_tcga_mutations' % c.lower()
                         for c in cancer_subtypes])
    parameters = {'cmd': 'getMutationData',
                  'gene_list': genes,
                  'genetic_profile_id': subtypes}

    r = requests.get(base_url, params=parameters)

    urlData = r.content
    error_message = 'Error: Problem when identifying'\
                    'a cancer study for the request.\n'
    if urlData == error_message:
        df = pd.read_table(io.StringIO(urlData.decode('utf-8')))
    else:
        df = pd.read_table(io.StringIO(urlData.decode('utf-8')), header=1)
        df = df[['gene_symbol', 'case_id',
                 'mutation_type', 'genetic_profile_id']]
        df = df[~((df.gene_symbol == 'Mutations') | (
            df.gene_symbol == 'gene_symbol'))]
        df = df.dropna()
    return df


def get_case_lists(cancer_subtypes):
    """ retrieve number of samples with mutations (missense, nonsense etc)
    and copy number abberations in the given list of cancer subtypes
    """
    base_url = 'http://www.cbioportal.org/webservice.do'

    subtypes = ' '.join(['%s_tcga' % c.lower()
                         for c in cancer_subtypes])
    parameters = {'cmd': 'getCaseLists',
                  'cancer_study_id': subtypes}
    r = requests.get(base_url, params=parameters)
    urlData = r.content
    df = pd.read_table(io.StringIO(urlData.decode('utf-8')))
    try:
        ng1 = df.case_list_description.iloc[2]
        result = re.search('samples \((.*) samples', ng1)
        num_mut_cases = int(result.group(1))
    except AttributeError:
        num_mut_cases = np.nan
    try:
        ng2 = df.case_list_description.iloc[3]
        result = re.search('CNA data \((.*) samples', ng2)
        num_cna_cases = int(result.group(1))
    except AttributeError:
        num_cna_cases = np.nan
    return num_mut_cases, num_cna_cases


def get_genetic_profile(cancer_subtypes):
    base_url = 'http://www.cbioportal.org/webservice.do'
    subtypes = ' '.join(['%s_tcga' % c.lower()
                         for c in cancer_subtypes])
    parameters = {'cmd': 'getGeneticProfiles',
                  'cancer_study_id': subtypes}
    r = requests.get(base_url, params=parameters)
    urlData = r.content
    df = pd.read_table(io.StringIO(urlData.decode('utf-8')))
    return df


def get_profile_data(gene_list, cancer_subtype):
    base_url = 'http://www.cbioportal.org/webservice.do'

    gene_list = sorted(gene_list)
    genes = ' '.join(gene_list)
    subtype = cancer_subtype.lower()

    parameters = {'cmd': 'getProfileData',
                  'case_set_id': "%s_tcga_cna" % subtype,
                  'genetic_profile_id': "%s_tcga_gistic" % subtype,
                  'gene_list': genes}

    r = requests.get(base_url, params=parameters)
    urlData = r.content
    error_message = 'Error: Problem when identifying a'\
                    'cancer study for the request.\n'
    if urlData == error_message:
        df = pd.read_table(io.StringIO(urlData.decode('utf-8')))
    else:
        columns = urlData.split('\n')
        columns2 = [c.split('\t')[2:] for c in columns[2:-1]]
        cols = ['Case_ID'] + gene_list
        ar = np.array(columns2).T
        df = pd.DataFrame(ar, columns=cols)
    return df


def compute_frequency_table(dark_kinases, subtypes,
                            mut_type='mut', cna_type='amplification',
                            subset=None):

    dfn = pd.read_csv('../data/dark_kinases/tcga_num_cases.csv')
    df_mut = pd.read_csv('../data/dark_kinases/dark_kinases_mutation_data.csv')
    if subset:
        df_mut = df_mut[df_mut.mutation_type.isin(subset)]

    total_mut = []
    for c in subtypes:
        kin_mut = []
        if mut_type == 'mut':
            total = dfn[dfn.Abbreviation == c].sequenced_samples.values[0]
            for kin in dark_kinases:
                dft = df_mut[(df_mut.gene_symbol == kin) & (
                    df_mut.genetic_profile_id ==
                    '%s_tcga_mutations' % c.lower())]
                if dft.empty:
                    mutated_samples = 0
                else:
                    mutated_samples = len(dft)
                try:
                    freq = 100 * mutated_samples/float(total)
                except ZeroDivisionError:
                    freq = 0
                kin_mut.append(freq)
        elif mut_type == 'cna':
            total = dfn[dfn.Abbreviation == c].samples_with_CNA.values[0]
            for kin in dark_kinases:
                try:
                    dfo = pd.read_csv('../data/dark_kinases/dark_kinases_'
                                      'profile_data_%s.csv' % c)
                    if cna_type == 'amplification':
                        nums = [s for s in dfo[kin].tolist() if s == 2]
                    elif cna_type == 'deletion':
                        nums = [s for s in dfo[kin].tolist() if s == -2]
                    mutated_samples = len(nums)
                except IOError:
                    mutated_samples = 0
                try:
                    freq = 100 * mutated_samples/float(total)
                except ZeroDivisionError:
                    freq = 0
                kin_mut.append(freq)
        total_mut.append(kin_mut)
    df = pd.DataFrame(np.array(total_mut),
                      columns=dark_kinases,
                      index=subtypes)
    return df
