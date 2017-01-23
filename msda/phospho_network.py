import pandas as pd
import networkx as nx
import requests
import numpy as np
import re
import subprocess
import mapping

file = ('/Users/kartik/Dropbox (HMS-LSP)/BrCaLines_profiling/'
        'RAWDATA/massspec/run2015/ReplicateA_pSTY_Summary_031315.xlsx')

df_ptm = pd.read_table('resources/Regulatory_sites_appended.tsv',
                       error_bad_lines=False)
df_kinase = pd.read_csv('resources/kinase_substrate_dataset_2016.csv')
# df_networkin = pd.read_table('resources/networkin_human_predictions.tsv')
df_networkin = pd.read_csv('resources/networkin_human'
                           '_predictions_appended.csv')


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
                            # 'siteIDstr': 'Site_Position',
                            'geneSymbol': 'Gene_Symbol',
                            'gene_symbol': 'Gene_Symbol',
                            'Gene Symbol': 'Gene_Symbol',
                            'motifPeptideStr': 'Motif',
                            'Localization score': 'Max_Score',
                            'Max Score': 'Max_Score'})
    return df


def split_sites(df, diff=None):
    """ refortm input dataframe by splitting phosphosite identifier (name_pos)
    to seperate columns for amino acid and position, retrieve necessary
     metadata for subsequent annotation with upstream kinases
    
    Parameter
    ---------
    df: pandas dataframe
        phosphoproteomics dataset with rows as features and samples as columns
    diff: str
       name of sample for which fold change is to be calculates

    Return
    ------
    df_clean: pandas dataframe
         dataframe that contains only metadata for all phosphosites
    """
    df = rename_columns(df)
    uids, names, motifs, sites, mx, type, fc = [], [], [], [], [], [], []
    for index in range(len(df)):
        motif = df.Motif.iloc[index]
        motif_list = motif.split(';')
        site = str(df.Site_Position.iloc[index])
        site_list = site.split(';')
        uids += [df.Uniprot_Id.iloc[index]] * len(motif_list)
        names += [df.Gene_Symbol.iloc[index]] * len(motif_list)
        try:
            type += [df.Origin.iloc[index].split('_')[1][0]] * len(motif_list)
        except IndexError:
            type += [df.Origin.iloc[index][0]] * len(motif_list)
        motifs += motif_list
        sites += site_list
        mx_score = str(df['Max_Score'].iloc[index]).split(';')
        mx += mx_score
        if diff is not None:
            fc += [df[diff].iloc[index]] * len(motif_list)
    try:
        uids = [id.split('|')[1] for id in uids]
    except IndexError:
        pass
    sites = ['%s%s' % (m[6], s) for m, s in zip(motifs, sites)]
    if diff is None:
        df_clean = pd.DataFrame(zip(uids, names, motifs, sites, mx, type),
                                columns=('Uniprot_Id', 'Gene_Symbol',
                                         'Motif', 'Site', 'score', 'type'))
    else:
        df_clean = pd.DataFrame(zip(uids, names, motifs, sites, mx, type, fc),
                                columns=('Uniprot_Id', 'Gene_Symbol',
                                         'Motif', 'Site', 'score', 'type', 'fc'))

    return df_clean


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
            print substrate
    df2 = pd.DataFrame(zip(ids, pos, aa))
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
    subprocess.call(['resources/NetworKIN_release3.0/NetworKIN.py',
                     '-n', 'resources/NetPhorest/netphorest',
                     '-b', 'resources/blast-2.2.17/bin/blastall', '9606',
                     fasfile, psitefile], stdout=f)
    return


def get_networkin_kinases(motif, df_nt):
    """ create list of kinase uniprot ids obtained 
    from Networkin prediction
    
    Parameter
    --------
    motif: str
        motif sequence of phosphosites
    df_nt: pandas dataframe
        dataframe of results from Netowkrin predcitons run locally

    Return
    ------
    kinase_uids: list
        list of kinase uniprto ids
    """
    motifp = motif[:5] + motif[5].lower() + motif[6:]

    precomputed_kinases = df_networkin[df_networkin.sequence == motifp][
        'string_path'].values.tolist()
    precomputed_kinases = [s.split(',')[0] for s in precomputed_kinases]
    precomp_scores = df_networkin[df_networkin.sequence == motifp][
        'networkin_score'].values.tolist()
    computed_kinases = df_nt[df_nt['Peptide sequence window'] == motifp][
        'Kinase/Phosphatase/Phospho-binding domain STRING ID'].values.tolist()
    comp_scores = df_nt[df_nt['Peptide sequence window'] == motifp][
        'NetworKIN score'].values.tolist()

    kinase_enps = precomputed_kinases + computed_kinases
    kinase_uids = [mapping.ensp2uid(id) for id in kinase_enps]
    networkin_scores = precomp_scores + comp_scores
    return kinase_uids, networkin_scores


def get_kinases(motif, organism=None):
    """ Get upstream kinases for given motif from PSP dataset

    Parameter
    ---------
    motif: str
       sequence of psp motif
    organism: str
       organism from which annotations are to be extracted

    Return
    ------
    kinases: list
       list of names of kinases
    kinase_ids: list
        list of kinase uniprot ids
    kin_orgs: list
        list of corrsponding organism for kinases
         in the kinase_substrate pair
    sub_orgs: list
        list of corresponding organism for substrates 
        in the kinase_substrate pair
    
    """
    if organism:
        df_org = df_kinase[df_kinase.SUB_ORGANISM == organism]
    else:
        df_org = df_kinase
    df_org['MOTIF'] = [mtf[2:-2].upper()
                       for mtf in df_org['SITE_+/-7_AA'].tolist()]

    kinases = df_org.KINASE[df_org.MOTIF == motif].values.tolist()
    kinase_ids = df_org.KIN_ACC_ID[df_org.MOTIF == motif].values.tolist()
    kin_orgs = df_org.KIN_ORGANISM[df_org.MOTIF == motif].values.tolist()
    sub_orgs = df_org.SUB_ORGANISM[df_org.MOTIF == motif].values.tolist()

    return kinases, kinase_ids, kin_orgs, sub_orgs


def generate_kinase_table(df_input_kinase, df_nt):
    """ Aggregate information on kinase annotations for 
    phosphosites from PSP and Networkin 

    Parameter
    ---------
    df_input: pandas dataframe
      kinase annotatiosn from Phosphosite
    df_nt: pandas dataframe
      kinase annotations form Networkin predictions

    Return
    -----
    df_output: pandas dataframe
       Long table of phosphosite metadata and corresponding kinases
    """
    df_input_kinase = df_input_kinase.drop_duplicates()
    substrate, site, ptype, kinases, kinase_ids = [], [], [], [], []
    source, motifs, confidence = [], [], []
    for ind, motif in enumerate(df_input_kinase.Motif.tolist()):
        out = get_kinases(motif)
        # kinases += out[0]
        kinase_ids += out[1]
        source += ['PSP'] * len(out[1])
        confidence += [100] * len(out[1])
        nkins = get_networkin_kinases(motif, df_nt)
        kinase_ids += nkins[0]
        source += ['Networkin'] * len(nkins[0])
        confidence += nkins[1]
        substrate += [df_input_kinase.Gene_Symbol.iloc[ind]] * (
            len(out[1]) + len(nkins[0]))
        site += [df_input_kinase.Site.iloc[ind]] * (
            len(out[1]) + len(nkins[0]))
        ptype += [df_input_kinase.type.iloc[ind]] * (
            len(out[1]) + len(nkins[0]))
        motifs += [motif] * (len(out[1]) + len(nkins[0]))
    
    site2 = ["%s_%s" % (s, t) for s, t in zip(site, ptype)]
    print len(substrate), len(source), len(site2), len(motifs), len(kinases)
    kinases = [mapping.uid2gn(id) for id in kinase_ids]
    df_output = pd.DataFrame(zip(substrate, site2, motifs,
                                 kinases, kinase_ids, source, confidence),
                             columns=['Gene_Symbol', 'Site', 'Motif',
                                      'KINASE', 'KINASE_ID', 'source',
                                      'confidence'])
    return df_output


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



def split_sites2(df, diff=None):
    """ refortm input dataframe by splitting phosphosite identifier (name_pos)
    to seperate columns for amino acid and position, retrieve necessary
     metadata for subsequent annotation with upstream kinases
    
    Parameter
    ---------
    df: pandas dataframe
        phosphoproteomics dataset with rows as features and samples as columns
    diff: str
       name of sample for which fold change is to be calculates

    Return
    ------
    df_clean: pandas dataframe
         dataframe that contains only metadata for all phosphosites
    """
    df = rename_columns(df)
    names, motifs, sites, type, fc = [], [], [], [], []
    for index in range(len(df)):
        motif = df.Motif.iloc[index][1:-1]
        motif_list = motif.split(';')
        site = str(df.Site_Position.iloc[index])
        site_list = site.split(';')
        # uids += [df.Uniprot_Id.iloc[index]] * len(motif_list)
        names += [df.Gene_Symbol.iloc[index]] * len(motif_list)
        type += [df.Origin.iloc[index][0]] * len(motif_list)
        motifs += motif_list
        sites += site_list
        # mx_score = str(df['Max_Score'].iloc[index]).split(';')
        # mx += mx_score
        if diff is not None:
            fc += [df[diff].iloc[index]] * len(motif_list)
    # try:
    #     uids = [id.split('|')[1] for id in uids]
    # except IndexError:
    #     pass
    sites = ['%s%s' % (m[6], s) for m, s in zip(motifs, sites)]
    if diff is None:
        df_clean = pd.DataFrame(zip(names, motifs, sites, type),
                                columns=('Gene_Symbol',
                                         'Motif', 'Site', 'type'))
    else:
        df_clean = pd.DataFrame(zip(names, motifs, sites, type, fc),
                                columns=('Gene_Symbol', 'Motif',
                                         'Site', 'type', 'fc'))

    return df_clean
