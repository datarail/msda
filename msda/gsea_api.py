import subprocess
import os


def get_gsea_enrichment(input_file, library, out_folder=None, out_file=None):
    path = "gsea_output/%s" % out_folder
    os.makedirs(path)
    lib = 'gsea_libraries/%s_v5.2.gmt' % library
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
