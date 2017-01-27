import subprocess
import os
import matplotlib.pyplot as plt
from pylab import barh
import numpy as np
import pandas as pd
from pylab import *
import seaborn as sns
import os


def get_gsea_enrichment(input_file, library, out_folder=None, out_file=None):
    path = "gsea_output/%s" % out_folder
    os.makedirs(path)
    lib = '%s.gmt' % library
    arg_list = ['java', '-cp', 'jars/gsea2-2.2.3.jar',
                '-Xmx1024m',  'xtools.gsea.GseaPreranked',
                '-param_file', 'gsea_parameters.txt',
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


def plot_nes(df, filter=False, top=5, outfile=None):
    df1 = df[df.NES > 0].iloc[:top, :]
    df2 = df[df.NES < 0].iloc[:top, :]
    df = pd.concat([df1, df2])
    if filter is True:
        df = df[(df['NOM p-val'] <= 0.05) & (df['FDR q-val'] <= 0.25)]
    plt.ion()
    f, ax = plt.subplots()
    ax = sns.barplot(x='NES', y='NAME', data=df, color='b')
    ax.set(ylabel="", xlabel="Normalzied Enrichment Score (NES)")
    pvals = ["%.3f" % p for p in df['NOM p-val'].tolist()]
    score_pos = [s+0.1 if s > 0 else s-0.3 for s in df['NES'].tolist()]
    ax.set_xlim([min(score_pos), max(score_pos)])
    for p, pval, sc in zip(ax.patches, pvals, score_pos):
        ya = p.get_y()
        ax.text(sc, ya+0.6, pval)
    if outfile:
        plt.savefig(outfile, dpi=600)
    else:
        plt.show()
    plt.clf()


def make_rnkfile(df, sample1, sample2):
    df2 = df.copy()
    df2[sample1] = df2[sample1].div(df2[sample2])
    df2[sample1] = df2[sample1].apply(np.log2)
    df2 = df2.sort([sample1])
    df2 = df2[['Gene_Symbol', sample1]]
    df2.Gene_Symbol = [g.upper() for g in df2.Gene_Symbol.tolist()]
    df2 = df2.replace([-np.inf, np.inf], np.nan)
    df2 = df2.dropna()
    return df2


