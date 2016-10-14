import pandas as pd
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
import networkx as nx

file = ('/Users/kartik/Dropbox (HMS-LSP)/BrCaLines_profiling/'
        'RAWDATA/massspec/run2015/ReplicateA_pSTY_Summary_031315.xlsx')

df_ptm = pd.read_table('resources/Regulatory_sites_appended.tsv',
                       error_bad_lines=False)
df_kinase = pd.read_csv('resources/kinase_substrate_dataset.csv')
df_networkin = pd.read_table('resources/networkin_human_predictions.tsv')

df1 = pd.read_excel(file, sheetname=0, header=1)
df2 = pd.read_excel(file, sheetname=1, header=1)
df3 = pd.read_excel(file, sheetname=2, header=1)
df4 = pd.read_excel(file, sheetname=3, header=1)


def rename_columns(df):
    df = df.rename(columns={'Protein Id': 'Protein_ID',
                            'proteinID': 'Protein_ID',
                            'Site Position': 'Site_Position',
                            # 'siteIDstr': 'Site_Position',
                            'geneSymbol': 'Gene_Symbol',
                            'gene_symbol': 'Gene_Symbol',
                            'motifPeptideStr': 'Motif'})
    return df


def split_sites(df):
    df = rename_columns(df)
    uids, names, motifs, sites = [], [], [], []
    for index in range(len(df)):
        motif = df.Motif.iloc[index]
        motif_list = motif.split(';')
        site = str(df.Site_Position.iloc[index])
        site_list = site.split(';')
        uids += [df.Protein_ID.iloc[index]] * len(motif_list)
        names += [df.Gene_Symbol.iloc[index]] * len(motif_list)
        motifs += motif_list
        sites += site_list
    uids = [id.split('|')[1] for id in uids]
    sites = ['%s%s' % (m[6], s) for m, s in zip(motifs, sites)]
    df_clean = pd.DataFrame(zip(uids, names, motifs, sites),
                            columns=('Protein_ID', 'Gene_Symbol',
                                     'Motif', 'Site'))
    return df_clean

df1 = split_sites(df1)
df2 = split_sites(df2)
df3 = split_sites(df3)
df4 = split_sites(df4)

df_input = pd.concat([df1, df2, df3, df4])

motifs = df_input.Motif.tolist()
uids = df_input.Protein_ID.tolist()

ptm_motifs = [m.upper()[1:-1] for m in df_ptm['SITE_+/-7_AA'].tolist()]
ptm_uids = df_ptm.ACC_ID.tolist()

kinase_motifs = [m.upper()[1:-1] for m in df_kinase['SITE_+/-7_AA'].tolist()]
kinase_uids = df_kinase.SUB_ACC_ID.tolist()

nkin_motifs = [m.upper() for m in df_networkin.sequence.tolist()]
kinase_motifs11 = [k[1:-1] for k in kinase_motifs]
pMS_motifs = [k[1:-1] for k in motifs]
df_input.Motif = pMS_motifs


venn2([set(motifs), set(ptm_motifs)],
      set_labels=('tnbc_phospho', 'downstream_annotations'))
plt.savefig('phosphoproteomics/dwnstrm_annotations_by_motif2.png')
plt.clf()

venn2([set(uids), set(ptm_uids)],
      set_labels=('tnbc_phospho', 'downstream_annotations'))
plt.savefig('phosphoproteomics/dwnstrm_annotaions_by_protein2.png')
plt.clf()

venn2([set(motifs), set(kinase_motifs)],
      set_labels=('tnbc_phospho', 'kinase_annotations'))
plt.savefig('phosphoproteomics/kinase_annotations_by_motif2.png')
plt.clf()

venn2([set(uids), set(kinase_uids)],
      set_labels=('tnbc_phospho', 'kinase_annotations'))
plt.savefig('phosphoproteomics/kinase_annotations_by_protein2.png')
plt.clf()

venn3([set(uids), set(kinase_uids), set(ptm_uids)],
      set_labels=('tnbc_phospho', 'kinase_annotations',
                  'downstream annotations'))
plt.savefig('phosphoproteomics/annotations_by_protein_venn3.png')
plt.clf()

venn3([set(motifs), set(kinase_motifs), set(ptm_motifs)],
      set_labels=('tnbc_phospho', 'kinase_annotations',
                  'downstream annotations'))
plt.savefig('phosphoproteomics/annotations_by_motif_venn3.png')
plt.clf()


def get_kinases(motif, organism=None):
    if organism:
        df_org = df_kinase[df_kinase.SUB_ORGANISM == organism]
    else:
        df_org = df_kinase
    df_org['MOTIF'] = [mtf[1:-1].upper()
                       for mtf in df_org['SITE_+/-7_AA'].tolist()]
    try:
        kinase = df_org.KINASE[df_org.MOTIF == motif].values[0]
        kinase_id = df_org.KIN_ACC_ID[df_org.MOTIF == motif].values[0]
        kin_org = df_org.KIN_ORGANISM[df_org.MOTIF == motif].values[0]
        sub_org = df_org.SUB_ORGANISM[df_org.MOTIF == motif].values[0]
    except IndexError:
        kinase = 'nan'
        kinase_id = 'nan'
        kin_org = 'nan'
        sub_org = 'nan'
    return kinase, kinase_id, kin_org, sub_org


def get_networkin_kinases(motif):
    motifp = motif[:6] + motif[6].lower() + motif[7:]
    df_motif = df_networkin[df_networkin.sequence == motifp[1:-1]]
    if not df_motif.empty:
        highest_score = df_motif.networkin_score.max()
        highest_scoring_kinase = df_motif.id[
            df_motif.networkin_score == highest_score].values[0]
    else:
        highest_scoring_kinase = 'nan'
    return highest_scoring_kinase


def generate_kinase_table(df_input):
    df_input_kinase = df_input[df_input.Motif.isin(all_motifs)]
    df_input_kinase = df_input_kinase.drop_duplicates()
    kinase_names, kinase_ids, orgs = [], [], []
    networkin_kinases = []
    for motif in df_input_kinase.Motif.tolist():
        out = get_kinases(motif)
        kinase_names.append(out[0])
        kinase_ids.append(out[1])
        orgs.append("%s_%s" % (out[2], out[3]))
        networkin_kinases.append(get_networkin_kinases(motif))
    df_output = df_input_kinase.copy()
    df_output['KINASE'] = kinase_names
    df_output['KINASE_ID'] = kinase_ids
    df_output['organisms'] = orgs
    df_output['networkin_kinase'] = networkin_kinases
    return df_output

df_out = generate_kinase_table(df_input)

subs = list(set(df_out.Gene_Symbol.tolist()))
kinases = list(set(df_out.KINASE.tolist()))
venn2([set(kinases), set(subs)], set_labels=('Kinases', 'Substrates'))
plt.savefig('phosphoproteomics/kinase_substrate.png')


def generate_network(df_output):
    G = nx.MultiDiGraph()
    for index in range(len(df_output)):
        kinase = df_output.KINASE.iloc[index]
        kinase_id = df_output.KINASE_ID.iloc[index]
        substrate = df_output.Gene_Symbol.iloc[index]
        sub_id = df_output.Protein_ID.iloc[index]
        site = df_output.Site.iloc[index]
        G.add_node(substrate, UP=sub_id)
        G.add_node(kinase, UP=kinase_id)
        G.add_edge(kinase, substrate, site=site)
    return G


def generate_ksea_library(kin_sub_table):
    df = pd.read_csv(kin_sub_table)
    psp_kinase = [k for k in df.KINASE.tolist() if k != float]
    networkin_kinases = [k for k in df.networkin_kinase.tolist()
                         if k != float]
    all_kinases = list(set(psp_kinase + networkin_kinases))

    gene_sets = []
    for kinase in all_kinases:
        df1 = df[df.KINASE == kinase]
        subs = df1.Gene_Symbol.tolist()
        sites = df1.Site.tolist()
        sub_sites = ['%s_%s' % (sub, site) for sub, site
                     in zip(subs, sites)]
        df2 = df[df.networkin_kinase == kinase]
        subs = df2.Gene_Symbol.tolist()
        sites = df2.Site.tolist()
        ss = ['%s_%s' % (sub, site) for sub, site in zip(subs, sites)]
        sub_sites = list(set(sub_sites + ss))
        if len(sub_sites) >= 25:
            gene_set = [kinase, ' '] + sub_sites
            gene_sets.append('\t'.join(gene_set))
    return gene_sets        
            
       
