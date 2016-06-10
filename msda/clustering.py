import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

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
    dendr = dendrogram(cell_line_clusters, labels=df_nrm.index)
    plt.tight_layout()
    plt.ylabel('Pearson correlation distance')
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
    plt.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.5)
    for sample, pc1, pc2 in zip(samples, X_pca[:,0], X_pca[:, 1]):
        plt.annotate(sample, xy=(pc1, pc2), textcoords='offset points')
    plt.savefig(plot_name)    
    plt.clf()

    
    
        
