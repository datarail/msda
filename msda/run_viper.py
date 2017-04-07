import subprocess
import os


def sum_duplicate_rows(df, identifier='Gene_Symbol',
                       filename='expr_temp.csv'):
    if identifier == 'Gene_Symbol':
        df[identifier] = [g.upper() for g in df[identifier].tolist()]
    df2 = df.sort_values(identifier).groupby(identifier).mean()
    df2.insert(0, identifier, df2.index.tolist())
    df2.to_csv(filename, index=False)
    return df2


def run_viper(df, type='ss', category=None, metafile=None, 
              test=None, ref=None, outfile='viper_out.csv',
              identifier='Gene_Symbol', regulon='brca-tf-regulon.rda'):
    df2 = sum_duplicate_rows(df, identifier)
    viper_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'viper')
    regulon_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'viper/regulons/%s' % regulon)
    if type == 'ss':
        subprocess.call(['Rscript', '--vanilla', os.path.join(viper_path, 'viper.R'),
                         'expr_temp.csv', metafile, outfile])
    elif type == 'ms':
        subprocess.call(['Rscript', '--vanilla', os.path.join(viper_path, 'msviper.R'),
                         'expr_temp.csv', metafile, category,
                         test, ref, outfile])
    elif type == 'fc':
        subprocess.call(['Rscript', '--vanilla', os.path.join(viper_path, 'fc_viper.R'),
                         'expr_temp.csv', regulon_path, outfile])

