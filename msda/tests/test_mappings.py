from msda import mapping

def test_get_uniprot_from_ensembl():
    id = mapping.get_uniprot_from_ensembl('ENSP00000266970')
    assert id == 'P24941'


def test_get_name_from_uniprot():
    gene_name = mapping.get_name_from_uniprot('P24941')
    assert gene_name == 'CDK2'


def test_get_entrez_from_name():
    entrez_id = mapping.get_entrez_from_name('CDK2')
    assert entrez_id == 1017


def test_get_name_from_entrez():
    gene_name = mapping.get_name_from_entrez(1017)
    assert gene_name == 'CDK2'


def test_get_uniprot_from_name():
    uid = mapping.get_uniprot_from_name('CDK2')
    assert uid == 'P24941'
    
