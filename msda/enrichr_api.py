import subprocess
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import time
import shutil


def prune_by_background(library, background):
    bc = open(background).readlines()
    bc_list = [pr.split('\n')[0] for pr in bc]
    terms = open('enrichr_libraries/%s.txt' % library).readlines()
    new_terms = []
    for term in terms:
        gene_list = term.split('\t')[2:-1]
        genes_in_background = [gene for gene in gene_list
                               if gene in bc_list]
        new_term = '\t'.join(term.split('\t')[0: 2] +
                             genes_in_background + ['\n'])
        new_terms.append(new_term)
    return new_terms


def get_enrichment(gene_list, library, background=None,
                   output_file='enrichr_output'):
    """ Run Enrichr given gene list and library
    
    Parameter
    ---------
    gene_list: list
        list of strings
    library: str
        Enrichr library without the file extension
    background: Bool
        True if library is to be pruned for platform background
    output_xml: str
    
    Return
    ------
    df: pandas dataframe
       dataframe containing enrichment results      
 
    """
    output_xml = "%s.xml" % output_file
    msda_path = os.path.dirname(os.path.abspath(__file__))
    enrichr_input_path = os.path.join(msda_path, 'enrichr_input.txt')
    
    with open(enrichr_input_path, 'w') as f:
        for gene in gene_list:
            f.write("%s\n" % gene)
    if background:
        lib_list = prune_by_background(library, background)
        lib_path = os.path.join(msda_path, 'enrichr_libraries_bc')
        with open('%s/%s.txt' % (lib_path, library), 'w') as f:
            for term in lib_list:
                f.write('%s\n' % term)        
    else:
        lib_path = os.path.join(msda_path, 'enrichr_libraries')
    lib = '%s/%s.txt' % (lib_path, library)    
    subprocess.call(['java', '-jar', os.path.join(msda_path, 'jars/enrichr.jar'),
                     enrichr_input_path, lib, output_xml])
    time.sleep(2)    
    if os.path.isfile('%s_%s.tsv' % (output_file, library)):
        df = pd.read_table('%s_%s.tsv' % (output_file, library))
        try:
            df = get_adjusted_pvals(df)
        except TypeError:
            pass
        # os.remove(lib)
    else:
        df = None
    return df


def get_adjusted_pvals(df):
    pvals = df['P-value'].tolist()
    adpvals = multipletests(pvals, 0.1, method='fdr_bh')
    df['Adjusted P-value'] = adpvals[1]
    return df
