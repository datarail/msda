def isnull_index(df, axis=0):
    ''' Given a pandas dataframe, return a dictionary that lists
    indeces of the missing values in each column '''

    assert axis < 2, "invalid axis dimension"

    if axis == 0:
        df = df.copy()
    elif axis == 1:
        df = df.transpose()
        
    isnull_idx = {}

    for col in df.columns:
        isnull_idx[col] = df[col][df[col].isnull()].index.tolist()

    return isnull_idx
