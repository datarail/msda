import subprocess


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
              identifier='Gene_Symbol'):
    df2 = sum_duplicate_rows(df, identifier)
    if type == 'ss':
        subprocess.call(['Rscript', '--vanilla', 'viper/viper.R',
                         'expr_temp.csv', metafile, outfile])
    elif type == 'ms':
        subprocess.call(['Rscript', '--vanilla', 'viper/msviper.R',
                         'expr_temp.csv', metafile, category,
                         test, ref, outfile])
    elif type == 'fc':
        subprocess.call(['Rscript', '--vanilla', 'viper/fc_viper.R',
                         'expr_temp.csv', outfile])

