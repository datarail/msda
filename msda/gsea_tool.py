import subprocess
import os
import matplotlib.pyplot as plt
from pylab import barh
import numpy as np
import pandas as pd
from pylab import *
import seaborn as sns


def get_gsea_enrichment(input_file, library, out_folder=None, out_file=None):
    path = out_folder
    os.makedirs(path)
    lib = '%s.gmt' % library
    msda_path = os.path.dirname(os.path.abspath(__file__))
    arg_list = ['java', '-cp', os.path.join(msda_path, 'jars/gsea2-2.2.3.jar'),
                '-Xmx1024m',  'xtools.gsea.GseaPreranked',
                '-param_file', os.path.join(msda_path, 'gsea_parameters.txt'),
                '-rnk', input_file,
                '-gmx', lib]
    if out_folder:
        arg_list += ['-out', path]
    if out_file:
        arg_list += ['-rpt_label', out_file]
    subprocess.call(arg_list)
    terminal_folder = os.listdir(path)[0]
    timestamp = terminal_folder.split('.')[2]
    terminal_path = "%s/%s" % (path, terminal_folder)

    file1 = '%s/gsea_report_for_na_neg_%s.xls' % (terminal_path, timestamp)
    file2 = '%s/gsea_report_for_na_pos_%s.xls' % (terminal_path, timestamp)
    try:
        df_pos = pd.read_table(file1)
        df_neg = pd.read_table(file2)
    except IOError:
        pass
    df = pd.concat([df_pos, df_neg])
    return df


def make_rnkfile(df, sample1, sample2, identifier):
    df2 = df.copy()
    df2[sample1] = df2[sample1].div(df2[sample2])
    df2[sample1] = df2[sample1].apply(np.log2)
    df2 = df2.sort([sample1])
    df2 = df2[[identifier, sample1]]
    df2[identifier] = [g.upper() for g in df2[identifier].tolist()]
    df2 = df2.replace([-np.inf, np.inf], np.nan)
    df2 = df2.dropna()
    return df2


def plot_nes(df, filter=False, top=10, outfile=None,
             pval_label='NOM p-val', fdr_label='FDR q-val',
             regulon='NAME',
             fontsize=12, show_pval=False):
    sns.set_style('dark')
    font = {'family': 'sans-serif',
            'color':  'black',
            'weight': 'normal',
            'size': fontsize,
            }
    sns.despine()
    if filter is True:
        df = df[(df[pval_label] <= 0.05) & (df[fdr_label] < 0.2)]
    df1 = df[df.NES > 0].iloc[:top, :]
    df2 = df[df.NES < 0].iloc[:top, :].sort_values(by='NES', ascending=False)
    ylen = np.max([len(df1), len(df2)])
    f, axes = plt.subplots(nrows=1, ncols=2)
    if not df1.empty:
        df1.NES.plot(ax=axes[1], kind='barh', width=0.3,
                     color='Green', alpha=0.5)
        for pos, name in zip(range(ylen), df1[regulon].tolist()):
            axes[1].text(0, pos-0.25, name, horizontalalignment='left',
                         fontdict=font)
    if not df2.empty:
        df2.NES.plot(ax=axes[0], kind='barh', width=0.3,
                     color='Red', alpha=0.5)
        for pos, name in zip(range(ylen), df2[regulon].tolist()):
            axes[0].text(0, pos+0.25, name, horizontalalignment='right',
                         fontdict=font)
    axes[0].set_yticks([])
    axes[1].set_yticks([])
    plt.subplots_adjust(wspace=0.075)
    axes[0].set_ylim([-0.5, ylen-0.5])
    axes[1].set_ylim([-0.5, ylen-0.5])

    if show_pval:
        pvals1 = ["%.2f" % p if p > 0 else "%.1f" % p
                  for p in df1[pval_label].tolist()]
        pvals2 = ["%.2f" % p if p > 0 else "%.1f" % p
                  for p in df2[pval_label].tolist()]
        for x_pos, y_pos, pval in zip(df1.NES.tolist(), range(ylen),
                                      pvals1):
            axes[1].text(x_pos+0.1, y_pos, pval)

        for x_pos, y_pos, pval in zip(df2.NES.tolist(), range(ylen)[::-1],
                                      pvals2):
            axes[1].text(x_pos-0.5, y_pos, pval)
    plt.gca().invert_yaxis()
    if not outfile:
        plt.show()
        plt.clf()
    else:
        plt.savefig(outfile)
        plt.clf()
