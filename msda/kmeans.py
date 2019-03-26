from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import matplotlib.cm as cm


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


#################################################################


def compare_silhoutte_scores(dfi, samples, range_n_clusters, cluster_dim='features'):
    """Compare silhoutte scores kmeans cluster numbers.
    Source code obtained and modified from :-
    http://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html

    Parameters
    ----------
    dfi : pandas dataframe
       input dataframe with features as rows and and samples as columns
    samples : list of str
       Names of samples
    range_n_clusters: list of int
      The list of cluster numbers for which the silhoutte score is to be computed.
    cluster_dim : Optional[str]
      Dimension along which data is to be clustered. Default is along features.
      To cluster samples, set cluster_dim='samples'

    Returns
    -------

    """
    df = dfi.fillna(0).copy()
    X = df[samples].values
    if cluster_dim == 'samples':
        X = X.T
    

    for n_clusters in range_n_clusters:
        # Create a subplot with 1 row and 2 columns
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(18, 7)

        # The 1st subplot is the silhouette plot
        # The silhouette coefficient can range from -1, 1 but in this example all
        # lie within [-0.1, 1]
        ax1.set_xlim([-0.1, 1])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

        # Initialize the clusterer with n_clusters value and a random generator
        # seed of 10 for reproducibility.
        clusterer = KMeans(n_clusters=n_clusters, random_state=10)
        cluster_labels = clusterer.fit_predict(X)

        # The silhouette_score gives the average value for all the samples.
        # This gives a perspective into the density and separation of the formed
        # clusters
        silhouette_avg = silhouette_score(X, cluster_labels)
        print("For n_clusters =", n_clusters,
              "The average silhouette_score is :", silhouette_avg)

        # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(X, cluster_labels)

        y_lower = 10
        for i in range(n_clusters):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = \
                sample_silhouette_values[cluster_labels == i]

            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.spectral(float(i) / n_clusters)
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                              0, ith_cluster_silhouette_values,
                              facecolor=color, edgecolor=color, alpha=0.7)

            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples

        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")

        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

        # 2nd Plot showing the actual clusters formed
        colors = cm.spectral(cluster_labels.astype(float) / n_clusters)
        ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                    c=colors, edgecolor='k')

        # Labeling the clusters
        centers = clusterer.cluster_centers_
        # Draw white circles at cluster centers
        ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
                    c="white", alpha=1, s=200, edgecolor='k')

        for i, c in enumerate(centers):
            ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
                        s=50, edgecolor='k')

        ax2.set_title("The visualization of the clustered data.")
        ax2.set_xlabel("Feature space for the 1st feature")
        ax2.set_ylabel("Feature space for the 2nd feature")

        plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                      "with n_clusters = %d" % n_clusters),
                     fontsize=14, fontweight='bold')

        plt.show()
