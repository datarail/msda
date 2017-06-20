from msda import mapping

def test_ensp2uid():
    id = mapping.ensp2uid('ENSP00000266970')
    assert id == 'P24941'


def test_uid2gn():
    gene_name = mapping.uid2gn('P24941')
    assert gene_name == 'CDK2'


def test_name2entrez():
    entrez_id = mapping.name2entrez('CDK2')
    assert entrez_id == 1017


def test_entrez2name():
    gene_name = mapping.entrez2name(1017)
    assert gene_name == 'CDK2'



def test_get_uniprot_id():
    uid = mapping.get_uniprot_id('CDK2')
    assert uid == 'P24941'
    
