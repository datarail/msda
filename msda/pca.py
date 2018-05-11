import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns
from matplotlib.patches import Ellipse
import pandas as pd
import numpy as np
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
from matplotlib.legend import Legend


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
                 color_col=None, color_dict=None, alpha_col=None,
                 size_col=None, size_scale=100,
                 sd_ellipse=False,
                 xmin=None, xmax=None,
                 ymin=None, ymax=None,
                 annotate_points=None):

    df_pca = dfpca.copy()
    fig, ax = plt.subplots()

    # Assign size of data points
    # --------------------------
    if size_col is None:
        size_list = 100  # plt.rcParams['lines.markersize'] ** 2
    else:
        smin = df_pca[size_col].min()
        smax = df_pca[size_col].max()
        size_list = [10 + (size_scale * (1 - (sl - smin)/(smax-smin)))
                     for sl in df_pca[size_col].tolist()]
    df_pca['size'] = size_list

    # Assign colors for data points
    # -----------------------------
    if color_col is None:
        df_pca['color'] = [(66/255, 134/255, 244/255)] * len(df_pca)
    elif color_dict is None:
        labels = df_pca[color_col].unique()
        colrs = sns.color_palette("husl", len(labels))
        color_dict = {l: c for l, c in
                      zip(labels, colrs)}
        df_pca['color'] = [color_dict[label] for label
                           in df_pca[color_col].tolist()]
    else:
        df_pca['color'] = [color_dict[label] for label
                           in df_pca[color_col].tolist()]

    # Assing alpha for data points
    # -----------------------------
    if alpha_col is None:
        df_pca['alpha'] = [0.8]*len(df_pca)
    else:
        amin = df_pca[alpha_col].min()
        amax = df_pca[alpha_col].max()
        alpha_list = [1 - (al - amin)/(amax-amin) for al
                      in df_pca[alpha_col].tolist()]
        df_pca['alpha'] = alpha_list

    df_pca['rgba'] = [rgb + (a,) for rgb, a in
                      zip(df_pca['color'], df_pca['alpha'])]

    x = df_pca['pca_1'].tolist()
    y = df_pca['pca_2'].tolist()
    ax.scatter(x, y,
               s=df_pca['size'].tolist(),
               facecolors=df_pca['rgba'].tolist(),
               edgecolors=df_pca['color'].tolist(),
               lw=1)
    ax.set_xlabel('PC 1 (%.2f%%)' % (100 * explained_variance[0]))
    ax.set_ylabel('PC 2 (%.2f%%)' % (100 * explained_variance[1]))
    if xmin is None:
        xmin = ax.get_xlim()[0]
    if xmax is None:
        xmax = ax.get_xlim()[1]
    if ymin is None:
        ymin = ax.get_ylim()[0]
    if ymax is None:
        ymax = ax.get_ylim()[1]
    ax.plot([0, 0], [ymin, ymax], '--k', alpha=0.5)
    ax.plot([xmin, xmax], [0, 0], '--k', alpha=0.5)
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((ymin, ymax))

    if sd_ellipse:
        cov = np.cov(x, y)
        lambda_, v = np.linalg.eig(cov)
        lambda_ = np.sqrt(lambda_)
        for j in [1, 2, 3]:
            ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                          width=lambda_[0] * j * 2,
                          height=lambda_[1] * j * 2,
                          angle=np.rad2deg(np.arccos(v[0, 0])))
            ell.set_facecolor('none')
            ell.set_edgecolor('black')
            ax.add_artist(ell)
    if color_col is not None:
        hue_recs = []
        for label in color_dict.keys():
            hue_recs.append(mpatches.Rectangle((0, 0), 1, 1,
                                               fc=color_dict[label]))
        legend1 = Legend(parent=ax, handles=hue_recs,
                         labels=color_dict.keys(), loc=(0.6, 0.02))
        ax.add_artist(legend1).get_frame().set_edgecolor("black")
    if annotate_points is not None:
        for i, ann in enumerate(df_pca[annotate_points].tolist()):
            ax.annotate(ann, (x[i], y[i]),
                        xytext=(-20, 20), textcoords='offset pixels',
                        fontsize=6)
    return df_pca


class HandlerEllipse(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = mpatches.Ellipse(xy=center, width=1.5*height + xdescent,
                             height=1.5*height + ydescent)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]
