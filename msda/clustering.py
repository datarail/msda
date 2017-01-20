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


def pca(df, meta_df, num_components=2, label=None, plot_prefix='pca_plot_'):
    df = df.copy()
    samples = meta_df.Sample.tolist()
    assert set(samples) < df.columns.tolist(), "sample names mismatched"
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
        plot_name = '%s%d_%d.png' % (plot_prefix, pcs[0]+1, pcs[1]+1)
        plot_pca(X_pca, explained_variance, samples, pcs, plot_name,
                 y, color_map)
    return X_pca, explained_variance


def plot_pca(X_pca, explained_variance, samples,
             pcs=[0, 1], plot_name='pca_plot_12.png', y=None, labels=None):
    Xs_pca = np.zeros([X_pca.shape[0], 2])
    Xs_pca[:, 0] = X_pca[:, pcs[0]]
    Xs_pca[:, 1] = X_pca[:, pcs[1]]
    if not labels:
        plt.scatter(Xs_pca[:, 0], Xs_pca[:, 1], alpha=0.5)  # s=20)
    elif labels:
        for lab, col in labels.iteritems():
            plt.scatter(Xs_pca[y == lab, 0], Xs_pca[y == lab, 1],
                        label=lab, color=col)
    # texts = []        
    for sample, pc1, pc2 in zip(samples, Xs_pca[:, 0], Xs_pca[:, 1]):
        plt.annotate(sample,
                     xy=(pc1, pc2), xytext=(-3, 3), size=5,
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
    plt.savefig(plot_name, dpi=600)
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


cell_line_response = {'A549': 'partial', 'Colo205': 'partial',
                      'DU145': 'partial', 'H1703': 'sensitive',
                      'HCT116': 'partial', 'HELA': 'partial',
                      'HT29': 'resistant', 'OVCAR8': 'resistant',
                      'SKBR3': 'resistant', 'SKOV3': 'resistant'}
