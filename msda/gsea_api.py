import subprocess
import os
import matplotlib.pyplot as plt
from pylab import barh
import numpy as np


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

    
def plot_nes(path, id):
    file1 = '%s/gsea_report_for_na_neg_%s.xls' % (path, id)
    file2 = '%s/gsea_report_for_na_pos_%s.xls' % (path, id)
    try:
        df1 = pd.read_table(file1)
        df2 = pd.read_table(file2)
    except IOError:
        pass
    names = df1.NAME.tolist()[:10]
    nes = df1.NES.tolist()[:10]
    names += df2.NAME.tolist()[:10]
    # names = '_'.join([n.split('_')[1:] for n in names])
    nes += df2.NES.tolist()[:10]
    barh(np.arange(len(nes)), nes)
    yticks(np.arange(len(nes)), tuple(names))
    lib = path.split('/')[-2].split('_')[-2]
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.show()
    plt.clf()

