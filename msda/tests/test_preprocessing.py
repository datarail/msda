import pandas as pd
from msda import preprocessing
from pandas.util.testing import assert_frame_equal


def test_quantile_normalize():
    df = pd.DataFrame({'C1': {'A': 5, 'B': 2, 'C': 3, 'D': 4},
                   'C2': {'A': 4, 'B': 1, 'C': 4, 'D': 2},
                   'C3': {'A': 3, 'B': 4, 'C': 6, 'D': 8}})
    dfq = preprocessing.quantile_normalize(df)
    dfq_ref = pd.DataFrame({'C1': {'A': 5.666667, 'B': 2.0,
                                   'C': 3.0, 'D': 4.666667},
                            'C2': {'A': 4.666667, 'B': 2,
                                   'C': 4.666667, 'D': 3.0},
                            'C3': {'A': 2.0, 'B': 3.0,
                                   'C': 4.666667, 'D': 5.666667}})
    
    assert_frame_equal(dfq, dfq_ref)
    

def test_filter_score():
    df = pd.DataFrame({'Max_Score': ['4;7', '15;2']},
                      index=['BRAF', 'RAF1'])
    dfs = preprocessing.filter_score(df)
    dfs_ref =  pd.DataFrame({'Max_Score': ['15;2']},
                      index=['RAF1'])
    assert_frame_equal(dfs, dfs_ref)


def test_get_primary_ids():
    secondary_id = 'A0A023LWJ7'
    primary_id = preprocessing.get_primary_ids(secondary_id)
    assert primary_id == 'Q0A255'
    
