import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
from matplotlib.gridspec import GridSpec
from matplotlib.legend import Legend
from matplotlib.legend_handler import HandlerPatch


def get_score(df, sample, signature):
    df = df.dropna(subset=[sample]).copy()
    fc = np.array(df[sample])
    genes = np.array(df.gene_name)
    fcp = fc[[x in signature for x in genes]]
    score = -np.mean(fcp)
    score = np.max((0.0012, score))
    return (np.log10(score) + 3)/3


def get_all_scores(df, samples, sig1, sig2):
    dfs = pd.DataFrame()
    for i, sample in enumerate(samples):
        if i % 10 == 0:
            print(i)
        score_dict = {}
        score_dict['sample'] = sample
        score_dict['sig1_score'] = get_score(df, sample, sig1)
        score_dict['sig2_score'] = get_score(df, sample, sig2)
        dfs = dfs.append(score_dict, ignore_index=True)
    return dfs


def plot_scatter(dfs, sig1, sig2, hue='agent', alpha='dose',
                 color_dict=None, marker='o', plot_legend=True,
                 ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    labels = dfs[hue].unique()
    if color_dict is None:
        colrs = sns.color_palette("husl", len(labels))
        color_dict = {}
        for i, label in enumerate(labels):
            color_dict[label] = colrs[i] + (1,)
    for ftr in dfs[hue].unique():
        dfis = dfs[dfs[hue] == ftr]
        x = dfis[sig1].tolist()
        y = dfis[sig2].tolist()
        rgba_colors = np.zeros((dfis.shape[0], 4))
        rgba_colors[:, :] = color_dict[ftr]
        edge_colors = rgba_colors.copy()
        if alpha is not None:
            amin = dfs[alpha].min()
            amax = dfs[alpha].max()
            rgba_colors[:, 3] = [1 - (al-amin)/amax for
                                 al in dfis[alpha].tolist()]
            sizes = [100 * al for
                     al in dfis[alpha].tolist()]
        ax.scatter(y, x, s=sizes, facecolor=rgba_colors,
                   edgecolor=edge_colors, lw=1, marker=marker)
    # Legends
    # -------
    if plot_legend:    
        hue_recs = []
        for label in color_dict.keys():
            hue_recs.append(mpatches.Rectangle((0, 0), 1, 1, fc=color_dict[label]))
        alpha_circles = []
        amin = dfs[alpha].min()
        amax = dfs[alpha].max()
        for ds in np.sort(dfs[alpha].unique()):
            alpha_val = 1 - (ds - amin)/amax
            gray_alpha = (128/255, 128/255, 128/255) + (alpha_val,)
            alpha_circles.append(
                mpatches.Circle((0, 0), radius=5, fc=gray_alpha))
        legend1 = Legend(parent=ax, handles=hue_recs,
                         labels=color_dict.keys(), loc=(0.6, 0.02), title='drug')
        legend2 = Legend(parent=ax, handles=alpha_circles,
                         labels=[1 * d for d in np.sort(dfs['dose'].unique())],
                         loc=4, title='dose',
                         handler_map={mpatches.Circle: HandlerEllipse()})
        ax.add_artist(legend2).get_frame().set_edgecolor("black")
        ax.add_artist(legend1).get_frame().set_edgecolor("black")


def get_kde(df, drug, dose, signature, grid=None, bandwidth=0.06):
    dg = df[(df.drug == drug) & (df.dose == dose)]
    # return dg[signature].tolist()
    if grid is None:
        grid = np.arange(0, 1.2, 0.02)
    x = np.array(dg[signature].tolist())
    kde = gaussian_kde(x, bandwidth / x.std(ddof=1))
    return kde.evaluate(grid)


def plot_kde(df, signature, cd=None, grid=None,
             ax=None, rotate=False, title=''):
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    if cd is None:
        cd = 'red'
    if grid is None:
        grid = np.arange(0, 1.2, 0.02)
    for drug in df.drug.unique():
        for dose in df.dose.unique():
            try:
                f_dist = get_kde(df, drug, dose, signature, grid)
                if rotate:
                    x, y = f_dist, grid
                else:
                    x, y = grid, f_dist
                ax.plot(x, y, c=cd[drug],
                        linewidth=np.max((1, dose)))
                if rotate:
                    ax.set_ylabel(title, fontweight='bold')
                    ax.set_ylim((0, df[signature].max() + 0.1))
                else:
                    ax.set_xlabel(title, fontweight='bold')
                    ax.set_xlim((0, df[signature].max() + 0.1))
                    ax.set_ylim((0, 3))
            except ValueError:
                print(drug, dose)


def cdk46_sl(df, samples, sig1, sig2):
    color_dict = {'Palbociclib': [0.97656, 0.19531, 0.19531, 1],
                  'Ribociclib': [0.8, 0.8, 0.5, 1],
                  'Abemaciclib': [0.13672, 0.23438, 0.99609, 1],
                  'Alvocidib': [0, 0, 0, 1]}
    gs = GridSpec(2, 2,
                  width_ratios=[1, 4],
                  height_ratios=[4, 1]
                  )
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax4 = plt.subplot(gs[3])

    dfs = get_all_scores(df, samples, sig1, sig2)
    dfs[['cell_line', 'drug', 'dose', 'time']] = dfs['sample'].apply(pd.Series)
    dfs['dose'] = [float(d) for d in dfs.dose.tolist()]
    dfs6 = dfs[(dfs.time == '6')]
    plot_scatter(dfs6, 'sig1_score', 'sig2_score',
                 hue='drug', alpha='dose', color_dict=color_dict, ax=ax2)
    dfs6s = dfs[(dfs.time == '6') & (dfs.dose != 1.0)]
    plot_kde(dfs6s, 'sig1_score', cd=color_dict, ax=ax1,
             rotate=True, title='G1-arrest score')
    plot_kde(dfs6s, 'sig2_score', cd=color_dict, ax=ax4, title='pan-CDK score')


def cdk46_dge(df, samples, sig1, sig2,
              rbnull=['BT549', 'PDX1258', 'PDXHCI002']):
    color_dict = {'Palbociclib': [0.97656, 0.19531, 0.19531, 1],
                  'Ribociclib': [0.8, 0.8, 0.5, 1],
                  'Abemaciclib': [0.13672, 0.23438, 0.99609, 1],
                  'Alvocidib': [0, 0, 0, 1]}
    gs = GridSpec(2, 2,
                  width_ratios=[1, 4],
                  height_ratios=[4, 1]
                  )
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax4 = plt.subplot(gs[3])

    dfs = get_all_scores(df, samples, sig1, sig2)
    dfs[['cell_line', 'drug', 'dose']] = dfs['sample'].str.split(
        '_', expand=True)
    dfs['dose'] = [round(float(d), 2) for d in dfs.dose.tolist()]
    dfs['dose'] = dfs['dose'].replace([3.17], 3.16)
    dfs1 = dfs[~dfs.cell_line.isin(rbnull)]  # [dfs.drug != 'Alvocidib']
    dfs2 = dfs[dfs.cell_line.isin(rbnull)]
    plot_scatter(dfs1, 'sig1_score', 'sig2_score',
                 hue='drug', alpha='dose', color_dict=color_dict,
                 ax=ax2, marker='o')
    plot_scatter(dfs2, 'sig1_score', 'sig2_score',
                 hue='drug', alpha='dose', color_dict=color_dict,
                 ax=ax2, marker='s',
                 plot_legend=False)
    dfs6s = dfs[(dfs.drug != 'Alvocidib') & (dfs.dose != 1.0)]
    plot_kde(dfs6s, 'sig1_score', cd=color_dict, ax=ax1,
             rotate=True, title='G1-arrest score')
    plot_kde(dfs6s, 'sig2_score', cd=color_dict, ax=ax4, title='pan-CDK score')


class HandlerEllipse(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = mpatches.Ellipse(xy=center, width=1.5*height + xdescent,
                             height=1.5*height + ydescent)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]
    
