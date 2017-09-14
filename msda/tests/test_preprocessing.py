import pandas as pd
import pytest
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


def test_remove_human_contaminants():
     df = pd.DataFrame({'Uniprot_Id': ['P1234', '2_HUMAN_contaminant',
                                       'Q4567', '4_HUMAN_contaminant'],
                        'C2': [4,  1, 4, 2],
                        'C3': [3, 5, 6, 8]},
                       index=[0, 1, 2, 3])
     
     dfs = preprocessing.remove_human_contaminants(df)
     dfs_ref = pd.DataFrame({'Uniprot_Id': ['P1234',  'Q4567'],
                        'C2': [4, 4],
                        'C3': [3, 6]},
                       index=[0, 2])
     
     assert_frame_equal(dfs, dfs_ref)


def test_remove_reverse_proteins():
     df = pd.DataFrame({'Uniprot_Id': ['P1234', '##P1234',
                                       'Q4567', '##Q4567'],
                        'C2': [4,  1, 4, 2],
                        'C3': [3, 5, 6, 8]},
                       index=[0, 1, 2, 3])
     
     dfs = preprocessing.remove_reverse_proteins(df)
     dfs_ref = pd.DataFrame({'Uniprot_Id': ['P1234',  'Q4567'],
                        'C2': [4, 4],
                        'C3': [3, 6]},
                       index=[0, 2])
     
     assert_frame_equal(dfs, dfs_ref)
    
    
def test_correct_gene_names():
     df = pd.DataFrame({'Gene_Symbol': ['CDK2', '2016-03-06 00:00:00'],
                        'Uniprot_Id': ['P24941',  'O60337'],
                        'C3': [3, 5]})
     dfg = preprocessing.correct_gene_names(df)
     dfg_ref = pd.DataFrame({'Gene_Symbol': ['CDK2', 'MARCH6'],
                        'Uniprot_Id': ['P24941',  'O60337'],
                        'C3': [3, 5]})
     assert_frame_equal(dfg, dfg_ref)


def test_correct_uniprot_identifiers():
     df = pd.DataFrame({'Uniprot_Id': ['sp|P24941-2|CDK2_HUMAN', 'sp|Q7Z469|HUMAN',
                                       'sp|A0A023GPJ5|HUMAN'],
                        'C2': [4,  1, 4]}, index=['A', 'B', 'C'])
     dfp = preprocessing.correct_uniprot_identifiers(df)
     dfp_ref =  pd.DataFrame({'Uniprot_Id': ['P24941-2','P13060'],
                              'Uniprot_Id_prefix': ['P24941','P13060'],
                        'C2': [4, 4]}, index=['A', 'C'])
     assert_frame_equal(dfp, dfp_ref)
                      
    


def test_get_primary_ids():
    secondary_id = 'A0A023LWJ7'
    primary_id = preprocessing.get_primary_ids(secondary_id)
    assert primary_id == 'Q0A255'


@pytest.mark.xfail(reason="Data file not checked in yet")
def test_get_activating_modifications():
    from msda import ptm_info
    mods = ptm_info.get_modifications('P15056')
    assert mods == ['S446-p', 'T599-p', 'S602-p', 'S579-p', 'S729-p']

    
@pytest.mark.xfail(reason="Data file not checked in yet")
def test_get_inhibitory_modifications():
    from msda import ptm_info
    mods = ptm_info.get_modifications('P06400', direction='inhibited')
    assert mods == ['S795-p', 'S780-p', 'S567-p']

    
    

    
