from mygene import MyGeneInfo


def get_uid(name):
    mg = MyGeneInfo()
    res = mg.query(name, scopes='symbol, alias',
                   fields='uniprot, symbol', species='human')
    symbol = []
    uid = []
    for hit in res['hits']:
        uid.append(hit['uniprot']['Swiss-Prot'])
        symbol.append(hit['symbol'])
    return zip(symbol, uid)

    


    
