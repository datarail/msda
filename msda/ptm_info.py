import pandas as pd


def get_regulatory_info(uid, site, mod):
    info = {}
    df = pd.read_csv('resources/ptm_mapping.csv')
    site_spec = "%s-%s" % (site, mod)
    truth_table = (df.Uniprot_ID == uid) & (df.Site.isin([site_spec]))
    functions = df.Function[truth_table].values[0]
    info['Function'] = functions.strip().split(';')
    processes = str(df.Process[truth_table].values[0])
    info['Process'] = processes.strip().split(';')
    info['PPI'] = df.PPI[truth_table].values[0]
    info['other'] = df.Other_interactions[truth_table].values[0]
    pmids = df.PMID[truth_table].values[0]
    info['evidence'] = pmids.strip().split(';')
    return info


def get_identifier(uid, site, mod, organism):
    df = pd.read_table('resources/phosphosite_dataset_appended.tsv')
    site_spec = "%s-%s" % (site, mod)
    id = df[(df.MOD_RSD == site_spec) &
            (df.ORGANISM == organism) &
            (df.ACC_ID == uid)].SITE_GRP_ID.values[0]
    return id


def generate_report(input_csv, organism):
    df_input = pd.read_csv(input_csv)
    df_output = df_input.iloc[:, 0:5].copy()
    df_output['PSP_idenitifier'] = ['NA']*len(df_output)
    df_output['Effect_on_function'] = ['NA']*len(df_output)
    df_output['Effect_on_process'] = ['NA']*len(df_output)
    df_output['Evidence'] = ['NA']*len(df_output)

    for ind in range(len(df_output)):
        uid = str(df_output.uniprotID.iloc[ind])
        aa = str(df_output['amino acid'].iloc[ind])
        site_num = df_input['Site Position'].ix[ind]
        mod = str(df_input['Mod_Type'].ix[ind])
        site = '%s%d' % (aa, site_num)
        try:
            id = get_identifier(uid, site, mod, organism)
            df_output['PSP_idenitifier'].iloc[ind] = id
            try:
                info = get_regulatory_info(uid, site, mod)
                df_output['Effect_on_function'].iloc[ind] = info['Function']
                df_output['Effect_on_process'].iloc[ind] = info['Process']
                df_output['Evidence'].iloc[ind] = info['evidence']
            except IndexError:
                pass
        except IndexError:
            pass
    return df_output    
 
