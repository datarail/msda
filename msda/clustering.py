from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from sklearn.lda import LDA


def normalize_min_max(df):
    df = df.copy()
    for feature in df.columns:
        min = df[feature].min()
        range = df[feature].max() - min
        df[feature] = (df[feature] - min) / range
    return df


def hierarchical_clustering(df, samples, plot_name='hc_plot.png'):
    df = df.copy()
    df = df.transpose()
    df = df.ix[samples]
    df_nrm = normalize_min_max(df)
    cell_line_clusters = linkage(
        pdist(df_nrm, metric='correlation'),
        method='complete')
    dendrogram(cell_line_clusters, labels=df_nrm.index)
    plt.ylabel('Pearson correlation distance')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(plot_name)
    plt.clf()


def pca(df, samples, plot_name='pca_plot.png'):
    df = df.copy()
    df = df.transpose()
    df = df.ix[samples]
    df_nrm = normalize_min_max(df)
    X = df_nrm.values
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)
    var1, var2 = pca.explained_variance_ratio_
    plt.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.5)
    for sample, pc1, pc2 in zip(samples, X_pca[:, 0], X_pca[:, 1]):
        plt.annotate(sample, xy=(pc1, pc2), textcoords='offset points')
    plt.xlabel('PCA1 (explained_variance = %.2f%%)' % (100*var1))
    plt.ylabel('PCA2 (explained_variance = %.2f%%)' % (100*var2))
    plt.savefig(plot_name)
    plt.clf()


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

cell_line_response = {'A549': 'partial', 'Colo205': 'partial',
                      'DU145': 'partial', 'H1703': 'sensitive',
                      'HCT116': 'partial', 'HELA': 'partial',
                      'HT29': 'resistant', 'OVCAR8': 'resistant',
                      'SKBR3': 'resistant', 'SKOV3': 'resistant'}
