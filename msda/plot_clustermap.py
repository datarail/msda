from matplotlib.colors import ListedColormap
import seaborn as sns
import pandas as pd
from itertools import chain
import os
import numpy as np


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


def get_metadata_colormap(dfm):
    df_cols = dfm.copy()
    color_dict = {}
    for feature in df_cols:
        labels = df_cols[feature].unique()
        if feature == 'dose':
            labels = np.sort([float(s) for s in labels])
            scheme = "YlOrRd"
        else:
            scheme = "husl"
        color_dict.update(construct_color_dict(labels, scheme))
        df_cols[feature] = df_cols[feature].map(color_dict)
    ocols = df_cols.columns.tolist()
    df_cols[''] = ['white'] * len(df_cols)
    ncols = [("%s " % s).split(' ') for s in ocols]
    nncols = list(chain.from_iterable(ncols))
    df_cols = df_cols[nncols]
    return df_cols, color_dict


def plot_clustermap(df_fc, dfm, yticklabels=False):
    dfc = df_fc.copy()
    dfc.columns = dfm.index.tolist()
    df_cols, color_dict = get_metadata_colormap(dfm)
    cmap = get_blue_yellow_cmap()
    cg = sns.clustermap(dfc.fillna(0), col_colors=df_cols,
                        cmap=cmap, vmin=-1.5, vmax=1.5,
                        yticklabels=yticklabels,
                        xticklabels=False)
    for label in color_dict.keys():
        cg.ax_col_dendrogram.bar(0, 0, color=color_dict[label],
                                 label=label, linewidth=0)
    cg.ax_col_dendrogram.legend(loc=(-0.6, -2), ncol=1)
