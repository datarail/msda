from matplotlib.colors import ListedColormap
import seaborn as sns
import pandas as pd
from itertools import chain
import os
import numpy as np
from msda import enrichr_api as ai
import matplotlib.pyplot as plt


def construct_color_dict(labels, scheme="husl"):
    colrs = sns.color_palette(scheme, len(labels)).as_hex()
    color_dict = {str(lab): col for lab, col
                  in zip(labels, colrs)}
    return color_dict


def get_blue_yellow_cmap():
    resource_path = os.path.join(
        os.path.dirname(
            os.path.abspath(__file__)
                       ), 'resources')
    df_hm = pd.read_csv(os.path.join(resource_path, "Heatmaps.tsv"),
                        header=None, sep='\t')
    hm = [tuple(x) for x in df_hm.values]
    cmap = ListedColormap(sns.color_palette(hm))
    return cmap


def get_metadata_colormap(dfm, features, color_dict=None):
    df_cols = dfm.copy()
    if color_dict is None:
        color_dict = {}
    for feature in features:
        labels = df_cols[feature].unique()
        keys_in_dict = list(color_dict.keys())
        labels = [l for l in labels if l not in keys_in_dict]
        if 'dose' in feature :
            labels = np.sort([float(s) for s in labels])
            scheme = "YlOrRd"
        else:
            scheme = "husl"
        color_dict.update(construct_color_dict(labels, scheme))
    for col in df_cols.columns.tolist():
        df_cols[col] = df_cols[col].map(color_dict)
    ocols = df_cols.columns.tolist()
    df_cols[''] = ['white'] * len(df_cols)
    ncols = [("%s " % s).split(' ') for s in ocols]
    nncols = list(chain.from_iterable(ncols))
    df_cols = df_cols[nncols]
    return df_cols, color_dict


def plot_clustermap(df_fc, df_cols, color_dict,
                    metric='euclidean', yticklabels=False):
    dfc = df_fc.copy()
    cmap = get_blue_yellow_cmap()
    cg = sns.clustermap(dfc.fillna(0), col_colors=df_cols,
                        cmap=cmap, vmin=-1.5, vmax=1.5,
                        metric=metric,
                        yticklabels=yticklabels,
                        xticklabels=False, figsize=(6, 7))
    for label in color_dict.keys():
        cg.ax_col_dendrogram.bar(0, 0, color=color_dict[label],
                                 label=label, linewidth=0)
    cg.ax_col_dendrogram.legend(loc=(-0.6, -2), ncol=1)
    return cg


def get_enriched_gene_list(df, cg, start, end):
    genes = df.index.tolist()
    reordered_ind = cg.dendrogram_row.reordered_ind
    reordered_genes = [genes[i] for i in reordered_ind]
    start_index = reordered_genes.index(start)
    end_index = reordered_genes.index(end) + 1
    gene_list = reordered_genes[start_index:end_index]
    return gene_list


def get_cmap_sig(gene_list):
    genes = [x for x in gene_list if x is not None]
    library = 'LINCS_L1000_Chem_Pert_down'
    dfe = ai.get_enrichment(genes, library)
    dfe = dfe[dfe['Adjusted P-value'] < 0.2]
    return dfe


def plot_cmap_sig(dfe):
    terms = []
    for t in dfe.Term.tolist():
        try:
            terms.append(t.split('-')[1])
        except IndexError:
            terms.append(t)
    pvals = [-np.log10(pv) for pv in dfe['Adjusted P-value']]
    dfp = pd.DataFrame(list(zip(terms, pvals)),
                       columns=['drug', '-log10(p-value)'])
    dfp2 = dfp[:10].copy()
    plt.barh(dfp2.index, dfp2['-log10(p-value)'], align='center')
    plt.yticks(dfp2.index, dfp2.drug)
    plt.gca().invert_yaxis()
    plt.subplots_adjust(left=0.2)
    plt.xlabel('False Discovery Rate')
    plt.title('L1000 drug_pert_down signature')
    return dfp
