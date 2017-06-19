import pandas as pd
import networkx as nx
import requests
import numpy as np
import re
import subprocess
import mapping
import os

resource_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'resources')
df_ptm = pd.read_table(os.path.join(resource_path,
                                    'Regulatory_sites_appended.tsv'),
                       error_bad_lines=False)
df_kinase = pd.read_csv(os.path.join(resource_path,
                                     'kinase_substrate_dataset_2016.csv'))
df_networkin = pd.read_csv(os.path.join(
    resource_path, 'networkin_human_predictions_appended.csv'))


def rename_columns(df):
    """ Rename columns so that they are standardized for further analysis
    Parameter
    ---------
    df: pandas dataframe

    Return
    ------
    df: pandas dataframe
    """

    df = df.rename(columns={'Protein Id': 'Uniprot_Id',
                            'proteinID': 'Uniprot_Id',
                            'Site Position': 'Site_Position',
                            'geneSymbol': 'Gene_Symbol',
                            'gene_symbol': 'Gene_Symbol',
                            'Gene Symbol': 'Gene_Symbol',
                            'motifPeptideStr': 'Motif',
                            'Localization score': 'Max_Score',
                            'Max Score': 'Max_Score'})
    return df


def get_annotated_subset(df_input):
    """ splits dataframe into subsets based on whether the 
    phosphosites have upstream kinases
    
    Parameter
    ---------
    df_input: pandas dataframe
       input dataframe with metadata for allphosphosites

    Return:
    df_annotated: pandas dataframe
       dataframe of phosphoite metadata for which annotation are 
       available on PSP and Networkin

    df_unnannotated: pandas dataframe
       dataframe pf phosphosite metdata for which annotations are not available
    """
    df_input.Motif = [m[1:-1] for m in df_input.Motif.tolist()]
    psp_motifs = [m.upper()[2:-2] for m  #
                  in df_kinase['SITE_+/-7_AA'].tolist()]
    nkin_motifs = [m.upper() for m
                   in df_networkin.sequence.tolist()]
    all_motifs = list(set(psp_motifs+nkin_motifs))
    df_annotated = df_input[df_input.Motif.isin(all_motifs)]
    df_unannotated = df_input[~df_input.Motif.isin(all_motifs)]
    return df_annotated, df_unannotated


def generate_network(df_output):
    """ generates kinase-substrate network using Networkx package
    """
    G = nx.MultiDiGraph()
    for index in range(len(df_output)):
        kinase = df_output.KINASE.iloc[index]
        kinase_id = df_output.KINASE_ID.iloc[index]
        substrate = df_output.Gene_Symbol.iloc[index]
        sub_id = df_output.Uniprot_Id.iloc[index]
        site = df_output.Site.iloc[index]
        G.add_node(substrate, UP=sub_id)
        G.add_node(kinase, UP=kinase_id)
        G.add_edge(kinase, substrate, site=site)
    return G


def generate_ksea_library(kin_sub_table, set_size=25):
    """ generate custom kinase set library with kinases as terms and 
    corresponding target  (phosphoites) list

    Parameter
    ---------
    kin_sub_table: csv file
       csv file that maps phosphosites to upstream kinases 
       from PSP and Networkin
    set_size: minimum size of kinase sets in the library

    Return
    ------
    gene_sets: list of strings
       each element in the list is a tab-separated entry of kinases 
       and downstream sites 
    """
    df = pd.read_csv(kin_sub_table)
    all_kinases = list(set([m.upper() for m in df.KINASE.tolist()]))
    gene_sets = []

    for kinase in all_kinases:
        df1 = df[df.KINASE == kinase]
        subs = [str(g).upper() for g in df1.Gene_Symbol.tolist()]
        sites = df1.Site.tolist()
        sub_sites = list(set(['%s_%s' % (sub, site) for sub, site
                              in zip(subs, sites)]))
        if len(sub_sites) >= set_size:
            gene_set = [kinase, ' '] + sub_sites
            gene_sets.append('\t'.join(gene_set))
    return gene_sets


def generate_substrate_fasta(df):
    """ gemerates fasta sequence files containing sequences of
    all proteins that contain phosphosites that do not have kinase
    annotations in PSP or Networkin. The outputs of the function
    will be used as input to run Networkin locally and predict kinases
   
    Parameter
    ---------
    df: pandas dataframe
       subset of phoproteomics data (metadata) that do
        not have kinase annotations

    Return
    ------
    substrate_fasta: list of strings
       each pair of elements in the list is a uniprot id (eg: '>P01345')
       followed by the sequence
    df2: pandas dataframe
       dataframe with uniprot id, amino acid and site of each phosphosite

    """

    substrate_fasta = []
    ids, aa, pos = [], [], []
    obsolete_entries = []
    for ind, substrate in enumerate(df.Uniprot_Id.tolist()):
        r = requests.get('http://www.uniprot.org/uniprot/%s.fasta' %
                         substrate)
        # substrate_fasta.append(r.text)
        seq_lines = r.text.split('\n')
        sequence = ''.join(seq_lines[1:])
        id_line = seq_lines[0]
        try:
            # id = re.search('>(.*)HUMAN', id_line).group(1) + 'HUMAN'
            id = re.search('>(?:sp|tr)\|(.*)\|', id_line).group(1)
            ids.append(id)
            # seq_lines[0] = id
            substrate_fasta.append(">%s\n%s\n" % (id, sequence))
            site = df.Site.iloc[ind]
            aa.append(site[0])
            pos.append(site[1:])
        except AttributeError:
            obsolete_entries.append(substrate)
    df2 = pd.DataFrame(zip(ids, pos, aa))
    if obsolete_entries:
        with open(os.path.join(resource_path, 'obsolete_entries.txt'), 'ab') as f:
            for s in list(set(obsolete_entries)):
                f.write("%s\n" % s)
    return substrate_fasta, df2


def create_rnk_file(df_input):
    """ creates file of fold change values to be used as input for GSEA

    Parameter
    ---------
    df_input: pandas dataframe
       input dataframe that contains phosphoste metadata and 
    fold change values (log2) of sample of interest

    Return
    ------
    df_rnk: pandas dataframe
       2-column dataframe of psp identifier (uniprotId_site) and 
    corresponding fold change value
    """
    fc = df_input.fc.tolist()
    gene = [str(g).upper() for g in df_input.Gene_Symbol.tolist()]
    site = df_input.Site.tolist()
    ty = [m[0].upper() for m in df_input.type.tolist()]

    id = ["%s_%s_%s" % (g, s, t) for g, s, t in zip(gene, site, ty)]
    df_rnk = pd.DataFrame(zip(id, fc), columns=('ps_id', 'fc'))
    df_rnk = df_rnk.sort(['fc'], ascending=True)
    return df_rnk


def run_networkin(fasfile, psitefile, outfile):
    """ run Networkin locally to generate kinase prediction
    Parameter:
    ----------
    fastafile: str
        path for sequence file (.fas) of all subtrates with unnanotated sites
    psitefile: str
        path for file (.res) with tab seperated values 
        of uniprot_id, amino acid and site
    outfile: str
        path for output file with Networkin predictions
    
    Return
    ------
    """
    f = open(outfile, 'wb')
    subprocess.call([os.path.join(resource_path, 'NetworKIN_release3.0/NetworKIN.py'),
                     '-n', os.path.join(resource_path, 'NetPhorest/netphorest'),
                     '-b', os.path.join(resource_path, 'blast-2.2.17/bin/blastall'), '9606',
                     fasfile, psitefile], stdout=f)
    return


def get_fc(df_input, samples, base_sample):
    """ compute fold change of samples relative to control
    Parameter
    ---------
    df_input: pandas dataframe
        phosphoproteomics dataset
    sample: list of str
       sample names
    base_sample:str
       sample that serves as control

    Return
    ------
    df: pandas dataframe
       samples normalized by base_sample
    """
    df = df_input[samples].div(df_input[base_sample], axis=0)
    df = df.apply(np.log2)
    return df


def series_split(df, col, type='str'):
    df[col] = [str(s) for s in df[col].tolist()]
    series = df[col].str.split(';').apply(pd.Series, 1).stack()
    series.index = series.index.droplevel(-1) # to line up with df's index
    series.name = col
    if type == 'float':
        series = series.astype(float)
    elif type == 'int':
        series = series.astype(float)
        series = series.astype(int)
    return series

def split_sites(df):

    df = rename_columns(df)
    origin = []
    for id in range(len(df)):
        if ';' in df.Motif.iloc[id]:
            origin.append('Composite')
        else:
            origin.append('Single')
    df['Origin'] = origin
    
    motif = series_split(df, 'Motif')
    max_score = series_split(df, 'Max_Score', type='float')
    sp = series_split(df, 'Site_Position', type='int')
    
    del df['Motif']
    del df['Max_Score']
    del df['Site_Position']
    dfe = pd.DataFrame(zip(motif, max_score, sp), index=sp.index.tolist(),
                       columns=[motif.name, max_score.name, sp.name])
    dfe['Site'] = ["%s%s" % (m[6], s) for m, s in zip(dfe.Motif.tolist(),
                                                      dfe.Site_Position.tolist())]
    df2 = df.join(dfe)
    #df2['Uniprot_Id'] = [u.split('|')[1] for u in df2.Uniprot_Id.tolist()]
    df2['Identifier'] = ["%s_%s%s_%s" % (gs, m[6], sp, o[0])
                         for gs, m, sp, o in zip(df2.Gene_Symbol.tolist(),
                                                  df2.Motif.tolist(),
                                                  df2.Site_Position.tolist(),
                                                  df2.Origin.tolist())]
    df2.index = df2.Identifier.tolist()
    return df2


def construct_table(df_nt, df_clean):
    col_rename = {'Target description': 'Gene_Symbol',
                  'Position': 'Site',
                  'Kinase/Phosphatase/Phospho-binding domain description': 'KINASE',
                  'NetworKIN score': 'confidence',
                  'Peptide sequence window': 'Motif'}
    df_nt = df_nt.rename(columns=col_rename)
    df_nt = df_nt[col_rename.values()]
    df_nt['Source'] = ['NetworKIN'] * len(df_nt)
    df_nt['Motif'] = [m.upper() for m in df_nt['Motif'].tolist()]
    df_nt.index = ["%s_%s" % (g, s) for g,s in zip(df_nt.Gene_Symbol.tolist(),
                                                   df_nt.Site.tolist())]
    df_kinase['SITE_+/-7_AA'] = [m[2:-2].upper() for m in df_kinase['SITE_+/-7_AA'].tolist()]
    df_clean.Motif = [m[1:-1] for m in df_clean.Motif.tolist()]
    df_psp = df_kinase[df_kinase['SITE_+/-7_AA'].isin(df_clean.Motif.tolist())]
    df_psp = df_psp.rename(columns={'SUBSTRATE': 'Gene_Symbol',
                                    'SUB_MOD_RSD': 'Site',
                                    'SITE_+/-7_AA': 'Motif'})
    df_psp['confidence'] = [100] * len(df_psp)
    df_psp = df_psp[col_rename.values()]
    df_psp['Source'] = ['PSP'] * len(df_psp)
    df_psp.index = ["%s_%s" % (g, s) for g,s in zip(df_psp.Gene_Symbol.tolist(),
                                                    df_psp.Site.tolist())]
    df_ks = pd.concat([df_nt, df_psp])
    return df_ks


def generate_kinase_annotations(df, path2data):
    """ Run networkin algorithm and lookup PSP to
        generate kinase-substrate library sets
    
    Parameter
    ---------
    df: pandas dataframe
       phosphoproteomics dataset 
    path2data: str
       local path to save annotation files generated
    
    
    Return
    ------
    df_out: pandas dataframe
       longtable of kinase annotatiosn from PSP and NetworKin
    """
    dfc = split_sites(df)
    subf, df_res = generate_substrate_fasta(dfc)

    fasfile = '%substrate.fas' % path2data
    resfile = '%spsite.res' % path2data
    outfile = '%sNetwork_predictions.txt' % path2data

    with open(fasfile, 'wb') as f:
        for line in subf:
            f.write(line)
    df_res.to_csv(resfile, index=False, sep='\t')

    # run networkin prediction
    run_networkin(fasfile, resfile, outfile)
    print "Running networkin algorithm. This may take 2+ hours"
    print "----------------------------------------------------"

    df_nt = pd.read_table(outfile)
    df_nt = df_nt[df_nt['NetworKIN score'] >= 2]
    df_out = construct_table(df_nt, df_clean)

    kin_table = '%skinase_substrate_table.csv' % path2data
    df_out.to_csv(kin_table, index=False)

    ksets = generate_ksea_library(kin_table, set_size=5)
    library = '%sksea_library' % path2data

    with open("%s.gmt" % library, 'wb') as f:
        for line in ksets:
            f.write("%s\n" % line)

    return df_out       
    
    

    
    
    
