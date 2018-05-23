from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import pandas as pd
import matplotlib.pyplot as plt
from msda import scatter


def compute_LDA(dfi, features, meta_col, num_components=2):
    """Ref http://python-for-multivariate-analysis.readthedocs.io/
       a_little_book_of_python_for_multivariate_analysis.html
       for useful examples of methods attributes.
    """
    # LDA
    dfs = dfi.dropna(subset=features)
    X = dfs[features].values
    y = dfs[meta_col].tolist()
    lda_cols = ['ld_%d' % ld for ld in range(1, num_components+1)]

    lda = LDA(n_components=num_components).fit(X, y)
    scalings = lda.scalings_
    df_loadings = pd.DataFrame(scalings[:, :num_components],
                               index=features, columns=lda_cols)

    explained_variance = lda.explained_variance_ratio_
    X_lda = lda.transform(X)

    df_lda = pd.DataFrame(X_lda)
    df_lda.columns = lda_cols
    df_lda.index = dfs.index.tolist()
    dfc = pd.concat([dfs, df_lda], axis=1)
    return dfc, explained_variance, df_loadings


def plot_scatter(dflda, explained_variance,
                 x_col='ld_1', y_col='ld_2',
                 color_col=None, color_dict=None, alpha_col=None,
                 size_col=None, size_scale=100,
                 sd_ellipse=False,
                 xmin=None, xmax=None,
                 ymin=None, ymax=None,
                 annotate_points=None):
    """ Scatter plot of LDA
    Parameters
    ----------
    dfpca : pandas dataframe
       Input dataframe that contains x and y linear discriminant
       columns used to plot as a scatter plot.
    explained_variance : list of floats
       the fractions of variance explained by each principal component.
    x_col : str
       column name corresponding to x-axis component
       of the scatter plot.
    y_col : str
       column name corresponding to y-axis compnent
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

    df_lda = dflda.copy()
    fig, ax = plt.subplots()
    cp1 = int(x_col[-1]) - 1
    cp2 = int(y_col[-1]) - 1
    xlabel = '%s (%.2f%%)' % (x_col.replace('_', '').upper(),
                              100 * explained_variance[cp1])
    ylabel = '%s (%.2f%%)' % (y_col.replace('_', '').upper(),
                              100 * explained_variance[cp2])
    df_lda = scatter.plot(df_lda, x_col=x_col, y_col=y_col,
                          color_col=color_col, color_dict=color_dict,
                          alpha_col=alpha_col,
                          size_col=size_col, size_scale=size_scale,
                          sd_ellipse=sd_ellipse,
                          xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                          annotate_points=annotate_points,
                          xlabel=xlabel, ylabel=ylabel, ax=ax)
    return df_lda


def plot_loading(df_loadings, component='ld_1', xlabel='Features'):
    cn = component.replace('_', ' ').upper()
    df_loadings = df_loadings.rename(columns={component: cn})
    dfl = df_loadings.sort_values(cn, ascending=False).copy()
    fig, ax = plt.subplots()
    dfl.plot(kind='bar', x=dfl.index, y=cn, color='blue', alpha=0.8, ax=ax)
    plt.subplots_adjust(bottom=0.25)
    xlim = plt.xlim()
    plt.plot(xlim, [0, 0], '--k', linewidth=1)
    plt.xlabel('Features', fontweight='bold')
    lead_features = dfl.index[0:5].tolist() + dfl.index[-5:].tolist()
    return lead_features
