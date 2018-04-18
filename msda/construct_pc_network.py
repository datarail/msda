import requests
import pandas as pd
import subprocess
import networkx as nx
import webbrowser
import subprocess
from py2cytoscape.data.cynetwork import CyNetwork
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.style import StyleUtil
import py2cytoscape.util.cytoscapejs as cyjs
import py2cytoscape.cytoscapejs as renderer



def get_pathsbetween(source_list, filter=False):
    """Pathway commons client that returns interactinos between
    genes in extended binary SIF format

    Parameters
    ----------
    source_list : list
        list of gene names
    filter : bool
        Set to TRUE to filter only those interactions that have the
        query genes as substrate and product

    Returns
    -------
    df : dataframe
        dataframe of interactions and supporting evidence
    """
    base_url = 'http://www.pathwaycommons.org/pc2/graph'
    source_str = ','.join(source_list)
    params = {'source': source_str,
              'kind': 'PATHSBETWEEN',
              'format': 'EXTENDED_BINARY_SIF'}
    r = requests.get(base_url, params=params)
    with open('temp_pc_interactions.txt', 'wb') as f:
        f.write(r.text)
    df = pd.read_table('temp_pc_interactions.txt')
    if filter:
        df = df[(df.PARTICIPANT_A.isin(source_list)) & (
            df.PARTICIPANT_B.isin(source_list))]
    return df


def make_network_plot(weights_file, network_file, figure_file, subsets=None):
    """Plots pathway commons network

    Parameters
    ----------
    weights_file : str
        path to file that contains 2-columns (genes and weights)
    network_file : str
        path to which output dataframe has to be saved
    figure_file : str
        path to which network figure is to be saved
    subsets : list
        list of gene names

    Return
    ------
    df2: dataframe
        3-colum dataframe

    """
    df = pd.read_csv(weights_file, sep='\t')
    source_list = df['0'].tolist()
    df = get_pathsbetween(source_list, filter=True)
    df.to_csv(network_file, index=False)
    df2 = df[df.columns[:3]]
    if subsets:
        df2 = df2[df2.INTERACTION_TYPE.isin(subsets)]
    df2.to_csv('temp_pc_file.txt', sep='\t', index=False)
    subprocess.call(['Rscript', 'netcontext.R', '-n', 'temp_pc_file.txt',
                     '-w',  weights_file, '-o', figure_file])
    return df2


def make_cytoscpe_plot(df_sif, weights_file):
    """Provide sif file and weights of nodes to launch and
    generate directed cytoscape network

    Parameters
    ----------
    df_sif : pandas dataframe
         sif file (can contain additional metadata on edge annotations)
    weights_file : str
        path to 2-column file of node weights
    """
    # Useful documentation
    # --------------------
    # http://nbviewer.jupyter.org/github/idekerlab/py2cytoscape/blob/develop/examples/New_wrapper_api_sample.ipynb
    cytoscape_path = '/Applications/Cytoscape_v3.4.0/Cytoscape.app'
    subprocess.call(['open', cytoscape_path])
    webbrowser.open('http://localhost:1234/v1')
    G = nx.from_pandas_dataframe(df_sif, 'PARTICIPANT_A',
                                 'PARTICIPANT_B', ['INTERACTION_TYPE'])
    df_weights = pd.read_csv(weights_file, sep='\t')
    df_weights.columns = ['Name', 'Weight']
    df_weights.index = df_weights.Name.tolist()

    cy = CyRestClient()
    cy.session.delete()
    cynet = cy.network.create_from_networkx(G)
    network_style = cy.style.create('GAL style')
    cy.layout.apply(network=cynet)
    cy.style.apply(network_style, network=cynet)
    cynet.update_node_table(df=df_weights, network_key_col='name')

    basic_settings = {

        # You can set default values as key-value pairs.

        'NODE_FILL_COLOR': '#6AACB8',
        'NODE_SIZE': 20,
        'NODE_BORDER_WIDTH': 0,
        'NODE_TRANSPARENCY': 120,
        'NODE_LABEL_COLOR': 'black',

        'EDGE_WIDTH': 2,
        'EDGE_TRANSPARENCY': 100,
        'EDGE_STROKE_UNSELECTED_PAINT': 'red',
        'EDGE_TARGET_ARROW_SHAPE': 'ARROW',
        'EDGE_TARGET_ARROW_UNSELECTED_PAINT': 'red',

        'NETWORK_BACKGROUND_PAINT': 'white'
    }

    network_style.update_defaults(basic_settings)
    network_style.create_passthrough_mapping(
        column='name', vp='NODE_LABEL', col_type='String')
    weights = cynet.get_node_column('Weight')
    color_gradients = StyleUtil.create_2_color_gradient(min=weights.min(),
                                                        max=weights.max(),
                                                        colors=('blue', 'red'))

    network_style.create_continuous_mapping(column='Weight',
                                            vp='NODE_FILL_COLOR',
                                            col_type='Double',
                                            points=color_gradients)

    edge_color_map = {'controls-state-change-of': 'orange',
                      'controls-phosphorylation-of': 'green',
                      'controls-transport-of': 'yellow',
                      'in-complex-with': 'blue',
                      'controls-expression-of': 'red'}

    edge_shape_map1 = {'controls-state-change-of': 'ARROW',
                       'controls-phosphorylation-of': 'ARROW',
                       'controls-transport-of': 'ARROW',
                       'in-complex-with': 'CIRCLE',
                       'controls-expression-of': 'ARROW'}

    edge_shape_map2 = {'controls-state-change-of': 'NONE',
                       'controls-phosphorylation-of': 'NONE',
                       'controls-transport-of': 'NONE',
                       'in-complex-with': 'CIRCLE',
                       'controls-expression-of': 'NONE'}

    network_style.create_discrete_mapping(column='INTERACTION_TYPE',
                                          col_type='String',
                                          vp='EDGE_STROKE_UNSELECTED_PAINT',
                                          mappings=edge_color_map)

    network_style.create_discrete_mapping(column='INTERACTION_TYPE',
                                          col_type='String',
                                          vp='EDGE_TARGET_ARROW_UNSELECTED_PAINT',
                                          mappings=edge_color_map)

    network_style.create_discrete_mapping(column='INTERACTION_TYPE',
                                          col_type='String',
                                          vp='EDGE_SOURCE_ARROW_UNSELECTED_PAINT',
                                          mappings=edge_color_map)

    network_style.create_discrete_mapping(column='INTERACTION_TYPE',
                                          col_type='String',
                                          vp='EDGE_TARGET_ARROW_SHAPE',
                                          mappings=edge_shape_map1)

    network_style.create_discrete_mapping(column='INTERACTION_TYPE',
                                          col_type='String',
                                          vp='EDGE_SOURCE_ARROW_SHAPE',
                                          mappings=edge_shape_map2)
