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
        df1 = pd.read_table(file1)
        df2 = pd.read_table(file2)
    except IOError:
        pass
    df = pd.concat([df1, df2])
    return df


def plot_nes(df, filter=False, outfile=None):
    if filter is True:
        df = df[(df['NOM p-val'] <= 0.05) & (df['FDR q-val'] <= 0.25)]
    # names = df.NAME.tolist()
    # names2 = [s.split('_') for s in names]
    # names3 = ['_'.join(s[1:]) for s in names2]
    # df.NAME = names3
    f, ax = plt.subplots()
    ax = sns.barplot(x='NES', y='NAME', data=df, color='b')
    ax.set(ylabel="", xlabel="Normalzied Enrichment Score (NES)")
    pvals = ["%.3f" % p for p in df['NOM p-val'].tolist()]
    score_pos = [s+0.1 if s > 0 else s-0.3 for s in df['NES'].tolist()]
    ax.set_xlim([min(score_pos), max(score_pos)])
    for p, pval, sc in zip(ax.patches, pvals, score_pos):
        # width = p.get_x()
        ya = p.get_y()
        ax.text(sc, ya+0.6, pval)
    if outfile:
        plt.savefig(outfile, dpi=600)
    else:
        plt.show()
    plt.clf()
    return df

    # nes = df.NES.tolist()
    # # names += df2.NAME.tolist()
    # # names = '_'.join([n.split('_')[1:] for n in names])
    # # nes += df2.NES.tolist()
    # # pos = np.arange(len(nes)) + 0.5
    # barh(pos, nes, align='center')
    # yticks(pos, tuple(names))
    # # lib = path.split('/')[-2].split('_')[-2]
    # plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    # xlabel('NES')
    # plt.savefig('nes1.png', dpi=600)
    # plt.clf()


def make_rnkfile(df, sample1, sample2):
    df2 = df.copy()
    df2[sample1] = df2[sample1].div(df2[sample2])
    df2[sample1] = df2[sample1].apply(np.log2)
    df2 = df2.sort([sample1])
    df2 = df2[['Gene_Symbol', sample1]]
    df2 = df2.replace([-np.inf, np.inf], np.nan)
    df2 = df2.dropna()
    return df2


