from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np


def cluster(dfi, samples, num_clusters=8, random_state=1):
    """Performs kmeans clustering on features in df

    Parameters
    ----------
    df : pandas dataframe
       dataframe of features (rna, protein, phospho) as rows
       and samples as columns
    samples : list of str
       Names of samples. Should match name in dataframe df.
    num_clusters : Optional[int]
       Number of clusters to form by kmeans
    random_state : Optional[int]
       Random seed for clustering should be defined for
       reproduabiity of results. Default is 1.

    Returns
    -------
    df : pandas dataframe
      input dataframe appended with cluster memebership column
    """
    df = dfi.fillna(0)
    X = df[samples].values
    kmeans = KMeans(n_clusters=num_clusters,
                    random_state=random_state).fit(X)
    cluster_number = kmeans.labels_
    df['kmeans_cluster_number'] = cluster_number
    return df


def plot(dfi, dfm, x_col, cluster_id_col='kmeans_cluster_number'):
    """Plots each cluster as trends across the samples.

    Parameters
    ----------
    dfi : pandas dataframe
      dataframe of features as rows and samples as columns. One column
      should contain info about the mapping of each feature to a cluster.
    dfm : pandas dataframe
      sample metadata table
    x_col : str
      Name of metadata column
    cluster_id_col : Optional[str]
      Name of column that contains cluster membership identities.

    Returns
    -------
    fig : matplotlib figure object

    """
    col_map = {c: m for c, m in
               zip(dfm['sample'], dfm[x_col])}
    df = dfi.rename(columns=col_map).copy()
    samples = list(col_map.values())
    clusters = df[cluster_id_col].unique()

    # Define plot properties
    # ----------------------
    grid_height = int(np.ceil(len(clusters) / 3))
    grid_dims = (grid_height, 3)
    fig = plt.figure(figsize=(35, 7 * grid_height), dpi=30)
    GridSpec(*grid_dims)

    # Loop over and plot each cluster as a subplot
    # --------------------------------------------
    for ai, cl in enumerate(clusters):
        dfc = df[df[cluster_id_col] == cl].copy()
        dfc = dfc[samples].mean(axis=1, level=0)
        ax_loc = np.unravel_index(ai, grid_dims)
        ax = plt.subplot2grid(grid_dims, ax_loc)
        dfc.T.plot(ax=ax, color='b',
                   xticks=range(dfc.shape[1]), legend=False)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.set_title(cl, fontsize=24, fontweight='bold')
    plt.subplots_adjust(hspace=0.35, bottom=0.1, top=0.95)
    return fig
