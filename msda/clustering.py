import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram

def normalize_min_max(df):
    df = df.copy()
    for feature in df.columns:
        min = df[feature].min()
        range = df[feature].max() - min
        df[feature] = (df[feature] - min) / range
    return df


def hierarchical_clustering(df, samples):
    df = df.copy()
    df = df.transpose()
    df = df.ix[samples]
    df_nrm = normalize_min_max(df)
    cell_line_clusters = linkage(
        pdist(df_nrm, metric='correlation'),
        method='complete')
    dendr = dendrogram(cell_line_clusters, labels=df_nrm.index)
    return dendr

    
        
