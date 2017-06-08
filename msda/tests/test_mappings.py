from msda import mapping

def test_ensp2uid():
    id = mapping.ensp2uid('ENSP00000266970')
    assert id == 'P24941'


def test_uid2gn():
    gene_name = mapping.uid2gn('P24941')
    assert gene_name == 'CDK2'



    
    
