import enrichr_client as ec
import requests
import logging


urllib3_logger = logging.getLogger('urllib3')
urllib3_logger.setLevel(logging.CRITICAL)

# enrichr_libraries = ['TRANSFAC_and_JASPAR_PWMs',
#                      'TargetScan_microRNA',
#                      'ENCODE_Histone_Modifications_2015',
#                      'KEGG_2016',
#                      'LINCS_L1000_Kinase_Perturbations_Down',
#                      'LINCS_L1000_Kinase_PErturbations_Up',
#                      'CORUM',
#                      'GO_Biological_Process_2015',
#                      'GO_Cellular_Component_2015',
#                      'GO_molecular_Function_2015',
#                      'LINCS_L1000_Chem_Pert_up',
#                      'LINCS_L1000_Chem_Pert_down',
#                      'LINCS_L1000_Ligand_Perturbations_up',
#                      'LINCS_L1000_Ligand_Perturbations_down',
#                      'OMIM_Disease',
#                      'Chromosome_Location',
#                      'ChEA_2015',
#                      'Pfam_InterPro_Domains']


enrichr_libraries = ['TRANSFAC_and_JASPAR_PWMs',
                     'KEGG_2016',
                     'GO_Biological_Process_2015',
                     'GO_Cellular_Component_2015',
                     'GO_molecular_Function_2015',
                     'ChEA_2015',
                     'ENCODE_TF_ChIP-seq_2015']

server = 'http://rest.ensembl.org'
ext = '/lookup/id'
headers = {"Content-Type": "application/json", "Accept": "application/json"}


def get_gene_name(ensembl_id):
    arg = '{"ids": ["%s"]}' % ensembl_id
    r = requests.post(server+ext, headers=headers, data=arg)
    try:
        r_dict = r.json()
        try:
            gene_name = r_dict[ensembl_id]['display_name']
        except TypeError:
            gene_name = 'trype_error'
        except KeyError:
            gene_name = 'NAN'
    except ValueError:
        gene_name = 'unknown'
    return gene_name


def make_gene_list(ensembl_list):
    ens_list = [e for e in ensembl_list if type(e) != float]
    gene_list = [get_gene_name(id) for id in ens_list]    
    return gene_list


def get_enichr_result(df):
    samples = df.columns.tolist()
    for sample in samples:
        gene_list = make_gene_list(df[sample].tolist())
        data = ec.post_gene_set(gene_list, sample)
        for library in enrichr_libraries:
            try:
                df_res = ec.get_enrichment_result(data, library)
            except AssertionError:
                print data['job_name'], library
            try:
                df_res.to_csv('../data/mouse/kelly_result2/%s_%s.csv' %
                              (sample, library),
                              index=False)
            except AttributeError:
                print data['job_name'], library


# for ens in ens3:
#     arg = '{"ids": ["%s"]}' % ens
#     r = requests.post(server+ext, headers=headers, data=arg)
#     try:
#         r_dict = r.json()
#         try:
#             gene_names3.append(r_dict[ens]['display_name'])
#         except TypeError:
#             gene_names3.append('unknown0')
#         except KeyError:
#             gene_names3.append('unknown_nan')
#     except ValueError:
#         gene_names3.append('unknown')




# data = ec.post_gene_set(gene_list, 'test')
# for library in enrichr_libraries:
#     df = ec.get_enrichment_result(data, library)
#     df.to_csv('../data/mouse/rna_up_minusF&P64<1_%s.csv' % library, index=False) 



