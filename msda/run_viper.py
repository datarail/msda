import subprocess


def sum_duplicate_rows(df, filename='expr_temp.csv'):
    df2 = df.sort_values('Gene_Symbol').groupby('Gene_Symbol').mean()
    df2.insert(0, 'Gene_Symbol', df2.index.tolist())
    df2.to_csv(filename, index=False)
    return df2


def run_viper(df, metafile, type='ss', category=None,
              test=None, ref=None, outfile='viper_out.csv'):
    d2 = sum_duplicate_rows(df)
    if type == 'ss':
        subprocess.call(['Rscript', '--vanilla', 'viper/viper.R',
                         'expr_temp.csv', metafile, outfile])
    elif type == 'ms':
        subprocess.call(['Rscript', '--vanilla', 'viper/msviper.R',
                         'expr_temp.csv', metafile, category,
                         test, ref, outfile])
