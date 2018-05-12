import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd
from msda import scatter


def compute_pca(df, dfm=None, num_components=2):
    if dfm is not None:
        samples = dfm['sample'].tolist()
    else:
        samples = df.columns.tolist()
    dfs = df[samples].dropna().transpose()
    X = dfs.values
    pca = PCA(n_components=num_components)
    X_pca = pca.fit_transform(X)
    df_pca = pd.DataFrame(X_pca)
    pca_cols = ['pca_%d' % pc for pc in range(1, num_components+1)]
    df_pca.columns = pca_cols
    df_pca.index = dfs.index.tolist()
    explained_variance = pca.explained_variance_ratio_
    if dfm is not None:
        dfm.index = dfm['sample'].tolist()
        df_pca = pd.concat([df_pca, dfm], axis=1)
        df_pca = df_pca.loc[:, ~df_pca.columns.duplicated()]
    return df_pca, explained_variance


def plot_scatter(dfpca, explained_variance,
                 x_col='pca_1', y_col='pca_2',
                 color_col=None, color_dict=None, alpha_col=None,
                 size_col=None, size_scale=100,
                 sd_ellipse=False,
                 xmin=None, xmax=None,
                 ymin=None, ymax=None,
                 annotate_points=None):

    df_pca = dfpca.copy()
    fig, ax = plt.subplots()
    xlabel = 'PC 1 (%.2f%%)' % (100 * explained_variance[0])
    ylabel = 'PC 2 (%.2f%%)' % (100 * explained_variance[1])
    df_pca = scatter.plot(df_pca, x_col='pca_1', y_col='pca_2',
                          xlabel=xlabel, ylabel=ylabel, ax=ax)
    return df_pca


