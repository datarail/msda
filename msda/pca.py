import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd
from msda import scatter


def compute_pca(df, dfm=None, num_components=2):
    """Performs PCA using sklearn's PCA function.

    Parameters
    ----------
    df : pandas dataframe
    dfm : Optional[pandas dataframe]
        metadata pertaining to samples (columns) in input dataframe.
    num_components : Optional[int]
       number of principal components to compute. Default is 2.
    Returns
    -------
    df_pca : pandas dataframe
       input dataframe dfm appended with the principal component vectors.
    explained_variance : list of floats
       the fractions of variance explained by each principal component.
    """
    if dfm is not None:
        samples = dfm['sample'].tolist()
    else:
        samples = df.columns.tolist()
    dfs = df[samples].dropna().transpose()
    X = dfs.values
    pca = PCA(n_components=num_components)
    X_pca = pca.fit_transform(X)
    df_pca = pd.DataFrame(X_pca)
    pca_cols = ['pc_%d' % pc for pc in range(1, num_components+1)]
    df_pca.columns = pca_cols
    df_pca.index = dfs.index.tolist()
    explained_variance = pca.explained_variance_ratio_
    if dfm is not None:
        dfm.index = dfm['sample'].tolist()
        df_pca = pd.concat([df_pca, dfm], axis=1)
        df_pca = df_pca.loc[:, ~df_pca.columns.duplicated()]
    return df_pca, explained_variance


def plot_scatter(dfpca, explained_variance,
                 x_col='pc_1', y_col='pc_2',
                 color_col=None, color_dict=None, alpha_col=None,
                 size_col=None, size_scale=100,
                 sd_ellipse=False,
                 xmin=None, xmax=None,
                 ymin=None, ymax=None,
                 annotate_points=None):
    """ Scatter plot of PCA

    Parameters
    ----------
    dfpca : pandas dataframe
       Input dataframe that contains x and y principal component
       columns used to plot as a scatter plot.
    explained_variance : list of floats
       the fractions of variance explained by each principal component.
    x_col : str
       column name corresponding to x-axis principal component
       of the scatter plot.
    y_col : str
       column name corresponding to y-axis principal compnent
       of the scatter plot.
    color_col : Optional[str]
       Name of metadata column that can be represented as different colors.
    color_dict : Optional[dict]
        key-value pairs where keys are unique labels in the color_col metadata
        column and values are rgb tuples
    alpha_col : Optional[str]
        Name of metadata column that can be represented on the data points plot
        as varying opacities.
    size_col : Optional[str]
        Name of metadata column that can be represented by size of data point.
    size_scale : Optional[float]
        Scaling factor to control the size of the data points.
    sd_ellipse : Optional[bool]
        Plots 1, 2 and 3 standard deviation envelopes.
    xmin : Optional[float]
        Minimum value of x-axis.
    xmax : Optional[float]
        Maximun value of x-axis.
    ymin : Optional[float]
        Minimum value of y-axis.
    ymax : Optional[float]
        Maximum value of y-axis
    annotate_points : Optional[str]
        Name of column containing metadata used to label the data points.

    Returns
    -------
    df_pca : pandas dataframe
       Input dataframe appended with color, opacity and size values
       of data points.
    """

    df_pca = dfpca.copy()
    fig, ax = plt.subplots()
    cp1 = int(x_col[-1]) - 1
    cp2 = int(y_col[-1]) - 1
    xlabel = '%s (%.2f%%)' % (x_col.replace('_', '').upper(),
                              100 * explained_variance[cp1])
    ylabel = '%s (%.2f%%)' % (y_col.replace('_', '').upper(),
                              100 * explained_variance[cp2])
    df_pca = scatter.plot(df_pca, x_col=x_col, y_col=y_col,
                          color_col=color_col, color_dict=color_dict,
                          alpha_col=alpha_col,
                          size_col=size_col, size_scale=size_scale,
                          sd_ellipse=sd_ellipse,
                          xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                          annotate_points=annotate_points,
                          xlabel=xlabel, ylabel=ylabel, ax=ax)
    return df_pca
