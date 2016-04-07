def pearsonr_pval(df, axis=0):
    ''' Given a pandas dataframe, the script computes pairwise correlation 
    and the associated p-values and stores the output in 2 dataframes '''

    from pandas import DataFrame
    from scipy.spatial.distance import pdist, squareform
    import scipy.stats as ss

    assert axis < 2, "invalid axis dimension."

    if axis == 1:
        df = df.transpose()
    elif axis == 0:
        df = df.copy()

    df_corr = DataFrame(squareform(
        pdist(df, lambda x, y: ss.pearsonr(x, y)[0])))

    df_pval = DataFrame(squareform(
        pdist(df, lambda x, y: ss.pearsonr(x, y)[1])))

    return df_corr, df_pval
