from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from sklearn.lda import LDA
import numpy as np
from itertools import combinations
from msda import adjustText
from collections import OrderedDict
import matplotlib.cm as cm
from collections import Counter
import seaborn as sns
import pandas as pd


def normalize_min_max(df):
    df = df.copy()
    for feature in df.columns:
        min = df[feature].min()
        range = df[feature].max() - min
        df[feature] = (df[feature] - min) / range
    return df


def hierarchical_clustering(df, samples, plot_name='hc_plot.png',
                            tl=True):
    df = df.copy()
    df = df.transpose()
    df = df.ix[samples]
    df_nrm = normalize_min_max(df)
    cell_line_clusters = linkage(
        pdist(df_nrm, metric='correlation'),
        method='complete')
    dendrogram(cell_line_clusters, labels=df_nrm.index,
               color_threshold=0)
    plt.ylabel('Pearson correlation distance')
    plt.xticks(rotation=90)
    if tl:
        plt.tight_layout()
    plt.savefig(plot_name, dpi=600)
    plt.clf()


def pca(df, meta_df, num_components=2, label=None, plot_prefix=None):
    df = df.copy()
    samples = number_duplicates(meta_df.Sample.tolist())
    samples2 = number_duplicates(df.columns.tolist())
    meta_df.Sample = samples
    df.columns = samples2
    df = df[samples]
    assert set(df.columns.tolist()) <= set(samples), "sample names mismatched"
    df = df.transpose()
    df = df.ix[samples]
    X = df.values
    if label:
        sample_map, color_map = generate_map(meta_df, label)
        y = np.array([sample_map[sample] for sample in df.index.tolist()])
    else:
        y = None
        color_map = None
    pca = PCA(n_components=num_components)
    X_pca = pca.fit_transform(X)
    explained_variance = pca.explained_variance_ratio_
    for pcs in list(combinations(range(num_components), r=2)):
        if plot_prefix:
            plot_name = '%s%d_%d.png' % (plot_prefix, pcs[0]+1, pcs[1]+1)
        else:
            plot_name = None
        plot_pca(X_pca, explained_variance, samples, pcs, plot_name,
                 y, color_map)
    return X_pca, explained_variance


def plot_pca(X_pca, explained_variance, samples,
             pcs=[0, 1], plot_name=None, y=None, labels=None):
    Xs_pca = np.zeros([X_pca.shape[0], 2])
    Xs_pca[:, 0] = X_pca[:, pcs[0]]
    Xs_pca[:, 1] = X_pca[:, pcs[1]]
    if not labels:
        plt.scatter(Xs_pca[:, 0], Xs_pca[:, 1], alpha=0.5, s=20)
    elif labels:
        for lab, col in labels.items():
            plt.scatter(Xs_pca[y == lab, 0], Xs_pca[y == lab, 1],
                        label=lab, color=col)
    # texts = []
    for sample, pc1, pc2 in zip(samples, Xs_pca[:, 0], Xs_pca[:, 1]):
        plt.annotate(sample, xy=(pc1, pc2), xytext=(-3, 3), size=5,
                     textcoords='offset points', ha='right', va='bottom')
#        texts.append(plt.text(pc1, pc2, sample,
#                              bbox={'pad': 0, 'alpha': 0}, size=7))
#    adjustText.adjust_text(Xs_pca[:,0], Xs_pca[:, 1], texts,
#                arrowprops=dict(arrowstyle="-", color='k', lw=0.5),
#                bbox={'pad':0, 'alpha':0}, size=7)

    plt.xlabel('PC %d (explained variance = %.2f%%)' %
               (pcs[0]+1, 100 * explained_variance[pcs[0]]))
    plt.ylabel('PC %d (explained variance = %.2f%%)' %
               (pcs[1]+1, 100 * explained_variance[pcs[1]]))
    plt.xticks([])
    plt.yticks([])
    plt.legend(fontsize=10, loc='lower left')
    # plt.tight_layout()
    if plot_name:
        plt.savefig(plot_name, dpi=600)
        plt.clf()
    else:
        plt.show()
        plt.clf()


def generate_map(meta_df, label):
    sample_map = OrderedDict()
    for sample in meta_df.Sample.tolist():
        sample_map[sample] = meta_df[label][
            meta_df.Sample == sample].values[0]
    sample_labels = list(set([types for types in sample_map.values()]))
    colors = cm.rainbow(np.linspace(0, 1, len(sample_labels)))
    label_map = {s: c for s, c in zip(sample_labels, colors)}
    return sample_map, label_map


def lda(df, samples, sample_labels, plot_name='lda_plot.png'):
    df = df.copy()
    df = df.transpose()
    df = df.ix[samples]
    df_nrm = normalize_min_max(df)
    X = df_nrm.values
    label_dict, y = encode_labels(sample_labels)
    ldas = LDA(n_components=2)
    X_lda = ldas.fit_transform(X, y)
    plot_scikit_lda(X_lda, y, label_dict, samples)


def encode_labels(sample_labels):
    enc = LabelEncoder()
    label_encoder = enc.fit(sample_labels)
    y = label_encoder.transform(sample_labels) + 1
    label_dict = {int: label for int, label in zip(y, sample_labels)}
    return label_dict, y


def plot_scikit_lda(X, y, label_dict, samples, plot_name='lda_plot.png'):

    ax = plt.subplot(111)
    for label, marker, color in zip(range(1, 4),
                                    ('^', 's', 'o'),
                                    ('blue', 'red', 'green')):

        plt.scatter(x=X[:, 0][y == label],
                    y=X[:, 1][y == label],  # flip the figure
                    marker=marker,
                    color=color,
                    alpha=0.5,
                    label=label_dict[label])

    for sample, ld1, ld2 in zip(samples, X[:, 0], X[:, 1]):
        plt.annotate(sample, xy=(ld1, ld2), textcoords='offset points')    
    plt.xlabel('LD1')
    plt.ylabel('LD2')

    leg = plt.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    # plt.title(title)

    # hide axis ticks
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                    labelbottom="on", left="off", right="off", labelleft="on")

    # remove axis spines
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)    

    plt.grid()
    plt.tight_layout
    plt.savefig(plot_name)
    plt.clf()


def diff_exp(df, samples, figname, top=25):
    df2 = df.copy()
    df2.index = df2.Gene_Symbol
    df2 = df2[samples]

    std = [np.std(df2.loc[protein].values)
           for protein in df2.index.values]
    df2['Standard_Deviation'] = std
    df2 = df2.sort('Standard_Deviation', ascending=False)
    df_nrm = normalize_min_max(df2[samples].transpose()).transpose()
    arr = df_nrm.ix[:top, samples].values
    plt.pcolor(arr, cmap='Blues', edgecolor='k')
    plt.gca().invert_yaxis()
    plt.yticks([i+0.5 for i in range(top)], df_nrm.index.values[:top])
    plt.xticks([i+0.5 for i in range(len(samples))], samples, rotation=90)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(figname)
    plt.clf()


def number_duplicates(samples):
    counts = Counter(samples)
    for s, num in counts.items():
        if num > 1:
            for suffix in range(1, num+1):
                samples[samples.index(s)] = s + str(suffix)
    return samples


def plot_clustermap(df, output_path, cmap=None, legend_label='',
                    z_score=None, xticklabels=False, yticklabels=True,
                    colors_dict=None, col_colors=None, row_colors=None):
    """Make clustermap figure

    Parameters
    ----------
    df
    df_meta
    output_path

    Returns
    -------
    cg
    """

    cg = sns.clustermap(df, col_colors=col_colors, row_colors=None,
                        cmap=cmap, z_score=z_score,
                        yticklabels=yticklabels, xticklabels=xticklabels)

    if colors_dict:
        for cat in colors_dict.keys():
            for label in colors_dict[cat]:
                cg.ax_col_dendrogram.bar(0, 0, color=colors_dict[cat][label],
                                         label=label, linewidth=0)
        cg.ax_col_dendrogram.legend(loc=(-0.7, -2), ncol=1)

    plt.subplots_adjust(top=1, bottom=0.02, left=0.3, right=0.8)
    fig = plt.gcf()
    fig.set_size_inches([10, 7.5])
    cg.cax.set_position((.025, .1, 0.025, .15))
    cg.cax.text(-0.3, -0.2, legend_label, fontsize=9)
    plt.savefig(output_path, dpi=300)
    return cg


def construct_categorical_pal(df_meta, category):
    categorical_pal = sns.color_palette("hls",
                                        len(df_meta[category].unique()))
    categorical_dict = dict(zip(df_meta[category].unique(), categorical_pal))
    return categorical_dict
