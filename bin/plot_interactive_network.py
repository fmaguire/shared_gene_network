#!/usr/bin/env python
from pyvis import network as net
import networkx as nx
import pickle
import argparse


def plot_interactive(networkx_pickle):
    network_width = "1080px"
    network_height = "1080px"

    nxg = nx.read_gpickle(networkx_pickle)

    # set attributes for visualisation from clearly labelled attributes
    clean_nx = nx.Graph()
    metadata_values = set()
    for node in nxg.nodes.data():
        if 'metadata' in node[1]:
            metadata_group = node[1]['metadata']
            metadata_values.update([metadata_group])
            clean_nx.add_node(node[0], title=f"{node[0].title()}: {metadata_group}", label=' ', group=metadata_group)
        else:
            clean_nx.add_node(node[0], title=node[0].title(), label=' ')


    max_edge_weight = max([x[2]['unique_genes'] for x in nxg.edges.data()])
    edge_padding = 5

    for edge in nxg.edges.data():
        edge_width = 2.5
        data = edge[2]
        if data['type'] == 'vf_only':
            clean_nx.add_edge(edge[0], edge[1], width=edge_width, color='orange', length=max_edge_weight - data['unique_genes'] + edge_padding, group=data['type'],
                          title=f"Unique Shared Genes: {data['unique_genes']}; AMR: {data['amr_set']}; VF: {data['vf_set']}")
        elif data['type'] == 'amr_only':
            clean_nx.add_edge(edge[0], edge[1], width=edge_width, color='green', length=max_edge_weight - data['unique_genes'] + edge_padding, group=data['type'],
                          title=f"Unique Shared Genes: {data['unique_genes']}; AMR: {data['amr_set']}; VF: {data['vf_set']}")
        elif data['type'] == 'amr_and_vf':
            clean_nx.add_edge(edge[0], edge[1], width=edge_width, color='purple', length=max_edge_weight - data['unique_genes'] + edge_padding, group=data['type'],
                          title=f"Unique Shared Genes: {data['unique_genes']}; AMR: {data['amr_set']}; VF: {data['vf_set']}")
        else:
            clean_nx.add_edge(edge[0], edge[1], width=edge_width, color='black', length=max_edge_weight - data['unique_genes'] + edge_padding, group=data['type'],
                          title=f"Unique Shared Genes: {data['unique_genes']}; AMR: {data['amr_set']}; VF: {data['vf_set']}")


    # add legends

    node_step = 70
    node_x = -int(network_width.replace('px', '')) / 2 + 50
    node_y = -int(network_height.replace('px', '')) /2 + 50

    if len(metadata_values) > 0:
        clean_nx.add_node('Node Legend', label="Node Legend", shape='box',
                          x=node_x, y=node_y, fixed=True, physics=False)
        for node in metadata_values:
            node_y += node_step
            clean_nx.add_node(node, label=node, group=node, x=node_x,
                              y=node_y, fixed=True, physics=False)


    edge_labels = {'green': 'AMR in shared genes',
                   'orange':'VF in shared genes',
                   'purple': 'AMR & VF in shared genes',
                   'black': 'Non-AMR/VF shared genes'}

    edge_x = int(network_width.replace('px', '')) / 2
    edge_y = -int(network_height.replace('px', '')) /2 + 50

    if len(metadata_values) > 0:
        edge_step = 25
    else:
        edge_step = 50
    clean_nx.add_node('Edge Legend', label="Edge Legend", shape='box',
                      x=edge_x, y=edge_y, fixed=True, physics=False)

    for colour, label in edge_labels.items():
        edge_y += edge_step
        clean_nx.add_node(label, label=label, shape='square',  color=colour,
                              x=edge_x, y=edge_y, fixed=True, physics=False)


    g = net.Network(height=network_height,
                    width=network_width,
                    heading="Shared Gene Network")

    g.set_options('''var options = {
        "edges": {
            "arrowStrikethrough": false,
            "color": {
                "inherit": true
            },
            "smooth": false
        },
        "interaction": {
            "hover": true,
            "navigationButtons": true,
            "tooltipDelay": 175
        }
    }
    ''')


    g.from_nx(clean_nx)
    g.write_html('network.html')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create interactive visualisation from pickled networkx graph')
    parser.add_argument('--networkx_pickle', help='Path to pickled networkx graph (generated in generate_graph.py)',
                        required=True)
    args = parser.parse_args()
    plot_interactive(args.networkx_pickle)

