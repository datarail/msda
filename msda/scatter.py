import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.legend import Legend
from matplotlib import colors
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def plot(df,
         x_col, y_col,
         color_col=None, color_dict=None,
         alpha_col=None, default_alpha=1,
         size_col=None, size_scale=100,
         sd_ellipse=False,
         xmin=None, xmax=None,
         ymin=None, ymax=None,
         annotate_points=None,
         xlabel=None, ylabel=None,
         ax=None, legend=True,
         default_color='b',
         rasterized=False):
    """Scatter plot based on input x and y columns in dataframe.
    Color, opacity and size of datapoints can also be set
    based on columns in the dataframe.

    Parameters
    ----------
    df : pandas dataframe
       Input dataframe that contains x and y columns
       used to plot as a scatter plot.
    x_col : str
       column name corresponding to x-axis of the scatter plot.
    y_col : str
       column name corresponding to y-axis of the scatter plot.
    color_col : Optional[str]
       Name of metadata column that can be represented as different colors.
    color_dict : Optional[dict]
        key-value pairs where keys are unique labels in the color_col metadata
        column and values are rgb tuples
    alpha_col : Optional[str]
        Name of metadata column that can be represented on the data points plot
        as varying opacities.
    default_alpha : Optional[int]
        If alpha column is not specified, alpha takes the value of defulat_alpha.
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
    xlabel : Optional[str]
        Label for x-axis.
    ylabel : Optional[str]
        Label for y-axis
    ax : Optional[subplot axes]
        Reference to subplot object.
    legend : Optional[bool]
        True if legend is to be shown. False otherwise. Default is True.
    default_color : Optional[str]
        Color of data points if color_dict is not provided. Default is blue.

    Returns
    -------
    dfs : pandas dataframe
        Input dataframe appended with values used for color, size and opacity.
    """

    dfs = df.copy()
    if ax is None:
        fig, ax = plt.subplots()

    # Assign size of data points
    # --------------------------
    if size_col is None:
        size_list = size_scale  # plt.rcParams['lines.markersize'] ** 2
    else:
        smin = dfs[size_col].min()
        smax = dfs[size_col].max()
        size_list = [10 + (size_scale * ((sl - smin)/(smax-smin)))
                     for sl in dfs[size_col].tolist()]
    dfs['size'] = size_list

    # Assign colors for data points
    # -----------------------------
    if color_col is None:
        rgb_color = colors.to_rgb(default_color)
        dfs['color'] = [rgb_color] * len(dfs)
    elif color_dict is None:
        labels = dfs[color_col].unique()
        colrs = sns.color_palette("husl", len(labels))
        color_dict = {l: c for l, c in
                      zip(labels, colrs)}
        dfs['color'] = [color_dict[label] for label
                        in dfs[color_col].tolist()]
    else:
        dfs['color'] = [color_dict[label] for label
                        in dfs[color_col].tolist()]

    # Assing alpha for data points
    # -----------------------------
    if alpha_col is None:
        dfs['alpha'] = [default_alpha] * len(dfs)
    else:
        amin = dfs[alpha_col].min()
        amax = dfs[alpha_col].max()
        alpha_list = [1 - (al - amin)/(amax-amin) for al
                      in dfs[alpha_col].tolist()]
        dfs['alpha'] = alpha_list

    dfs['rgba'] = [rgb + (a,) for rgb, a in
                   zip(dfs['color'], dfs['alpha'])]
    

    if 'markers' in dfs.columns.tolist():
        unique_markers = dfs['markers'].unique()
        for um in unique_markers:
            dfsm = dfs[dfs.markers == um].copy()          
    
            x = dfsm[x_col].tolist()
            y = dfsm[y_col].tolist()
            ax.scatter(x, y,
                       s=dfsm['size'].tolist(),
                       facecolors=dfsm['rgba'].tolist(),
                       edgecolors=dfsm['color'].tolist(),
                       marker=um,
                       lw=1, rasterized=rasterized)
            
    else:
        x = dfs[x_col].tolist()
        y = dfs[y_col].tolist()
        ax.scatter(x, y,
                   s=dfs['size'].tolist(),
                   facecolors=dfs['rgba'].tolist(),
                   edgecolors=dfs['color'].tolist(),
                   lw=1, rasterized=rasterized)

    if xlabel is None:
        xlabel = x_col
    if ylabel is None:
        ylabel = y_col
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')

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
    if legend & (color_col is not None):
        hue_recs = []
        for label in color_dict.keys():
            hue_recs.append(mpatches.Rectangle((0, 0), 1, 1,
                                               fc=color_dict[label]))
        legend1 = Legend(parent=ax, handles=hue_recs,
                         labels=color_dict.keys(), loc=(0.6, 0.02))
        ax.add_artist(legend1).get_frame().set_edgecolor("black")
    if annotate_points is not None:
        for i, ann in enumerate(dfs[annotate_points].tolist()):
            ax.annotate(ann, (x[i], y[i]),
                        xytext=(-20, 20), textcoords='offset pixels',
                        fontsize=6)
    return dfs
