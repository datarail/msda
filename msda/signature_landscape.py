import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde


def get_score(df, sample, signature):
    df = df.dropna(subset=[sample]).copy()
    fc = np.array(df[sample])
    genes = np.array(df.gene_name)
    fcp = fc[[x in signature for x in genes]]
    score = -np.mean(fcp)
    return score


def get_all_scores(df, samples, sig1, sig2):
    dfs = pd.DataFrame()
    for sample in samples:
        score_dict = {}
        score_dict['sample'] = sample
        score_dict['sig1_score'] = get_score(df, sample, sig1)
        score_dict['sig2_score'] = get_score(df, sample, sig2)
        dfs = dfs.append(score_dict, ignore_index=True)
    return dfs

    # color_dict = {'Palbociclib': [0.97656, 0.19531, 0.19531, 1],
    #              'Ribociclib': [0.8, 0.8, 0.5, 1],
    #              'Abemaciclib': [0.13672, 0.23438, 0.99609, 1]}


def plot_scatter(dfs, sig1, sig2, hue='cell_line', alpha='dose',
                 color_dict=None):
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
            amin = dfis[alpha].min()
            amax = dfis[alpha].max()
            rgba_colors[:, 3] = [1 - (al-amin)/amax for
                                 al in dfis[alpha].tolist()]
        plt.scatter(x, y, s=300, facecolor=rgba_colors,
                    edgecolor=edge_colors, lw=1)
    hue_recs = []
    for label in color_dict.keys():
        hue_recs.append(mpatches.Rectangle((0, 0), 1, 1, fc=color_dict[label]))
    alpha_circles = []
    amin = dfs[alpha].min()
    amax = dfs[alpha].max()
    for ds in np.sort(dfs[alpha].unique()):
        alpha_val = 1 - (ds - amin)/amax
        gray_alpha = (128/255, 128/255, 128/255) + (alpha_val,)
        print(gray_alpha)
        alpha_circles.append(
            mpatches.Circle((0, 0), fc=gray_alpha))
    ax.set_xlabel(sig1)
    ax.set_ylabel(sig2)
    legend1 = plt.legend(hue_recs, color_dict.keys(), loc=4)
    legend2 = plt.legend(alpha_circles, np.sort(dfs[alpha].unique()), loc=2)
    ax.add_artist(legend1)
    ax.add_artist(legend2)


def get_kde(df, sample, signature, grid=None):
    df = df.dropna(subset=[sample]).copy()
    fc = np.array(df[sample])
    genes = np.array(df.gene_name)
    fcp = fc[[x in signature for x in genes]]
    if grid is None:
        grid = np.arange(0, 1.2, 0.02)
    kde = gaussian_kde(fcp)
    return kde.evaluate(grid)
