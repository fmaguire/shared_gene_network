#!/usr/bin/env python

import pandas as pd
import pickle
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import itertools
import sys

def run(args):
    # parse RGI and index with simplified ORF name for lookup
    rgi_df  = pd.read_csv(args.rgi_output, sep='\t')
    rgi_df['ORF'] = rgi_df['ORF_ID'].str.split(' ').str.get(0)

    # similarly for VFDB results
    vfdb_df = pd.read_csv(args.vf_blast_output, sep='\t',
            names='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split())

    # finally grab the proteinortho output
    poff_df = pd.read_csv(args.proteinortho_output, sep='\t')
    poff_df = poff_df.rename(columns=lambda x: x.replace('.faa', ''))
    # remove any ORF that is only in one sample
    poff_df = poff_df[poff_df['# Species'] > 1]

    if args.seq_metadata:
        metadata_df = pd.read_csv(args.seq_metadata, sep='\t')
        metadata_df = metadata_df.set_index('name')
        metadata_annotate = metadata_df[args.metadata_col].to_dict()
    else:
        metadata_annotate = None

    # create graph with each plasmid as a node
    G = nx.Graph()
    for node in [x for x in poff_df.columns[3:]]:
        if metadata_annotate:
            if node in metadata_annotate:
                G.add_node(node, metadata = metadata_annotate[node])
            else:
                print(f"WARNING: {node} not present in 'name' column of {args.seq_metadata}", file=sys.stderr)
                G.add_node(node, metadata = "None")
        else:
            G.add_node(node)
    G.add_nodes_from([x for x in poff_df.columns[3:]])

    for index, ortho in poff_df.iterrows():
        # drop plasmids without that orthologue
        ortho = ortho[ortho != '*'].drop(['# Species', "Genes", "Alg.-Conn."]).to_dict()

        # if a single node edge has two # Species for some reason skip to next row
        if len(ortho) < 2:
            continue

        # extract all edges that orth implies
        all_edges = list(itertools.combinations(ortho.keys(), 2))

        all_orfs = list(itertools.chain(*[x.split(',') for x in ortho.values()]))

        # if any of the ORFs making up that orthalogue are an AMR gene then
        # mark the edge as containing an AMR genes (colour if any of the shared
        # genes are AMR
        amr_hits = set(rgi_df[rgi_df['ORF'].isin(all_orfs)]['Best_Hit_ARO'].values)
        vf_hits = set(vfdb_df[vfdb_df['qseqid'].isin(all_orfs)]['sseqid'].values)
        other_genes = set(all_orfs) - set(rgi_df[rgi_df['ORF'].isin(all_orfs)]['ORF']) - set(vfdb_df[vfdb_df['qseqid'].isin(all_orfs)]['qseqid'])

        # for each edge combination increase edge weight by one
        # i.e. number of unique genes shared between nodes/contigs/plasmids

        for edge in all_edges:
            # add edge if it doesn't already exist create it as null edge
            if not G.has_edge(*edge):
                G.add_edge(*edge, unique_genes=0, amr_set = set(),
                            vf_set = set(), other_genes = set())

            G.edges[edge]['unique_genes'] += 1

                # if not already true set amr status of edge
            if len(amr_hits) > 0:
                G.edges[edge]['amr_set'].update(amr_hits)

            if len(vf_hits) > 0:
                G.edges[edge]['vf_set'].update(vf_hits)

            if len(other_genes) > 0:
                G.edges[edge]['other_genes'].update(other_genes)

    # colour edges green if shared genes between plasmids include an AMR gene
    # dash the line if shared gene between plasmids include a VF gene
    all_edges = set([(u, v) for (u, v, d) in G.edges(data=True)])

    amr = set([(u, v) for (u, v, d) in G.edges(data=True) if len(d["amr_set"]) > 0])
    vf = set([(u, v) for (u, v, d) in G.edges(data=True) if len(d['vf_set']) > 0 ])
    both = list(amr.intersection(vf))
    amr_only = list(amr - vf)
    vf_only = list(vf - amr)
    neither = list(all_edges - vf - amr)

    # edges
    for edge in both:
        G.edges[edge]['type'] = 'amr_and_vf'
    for edge in amr_only:
        G.edges[edge]['type'] = 'amr_only'
    for edge in vf_only:
        G.edges[edge]['type'] = 'vf_only'
    for edge in neither:
        G.edges[edge]['type'] = 'neither'

    # change set to lists to allow serialising
    for edge in G.edges:
        G.edges[edge]['amr_set'] = list(G.edges[edge]['amr_set'])
        G.edges[edge]['vf_set'] = list(G.edges[edge]['vf_set'])
        G.edges[edge]['other_genes'] = list(G.edges[edge]['other_genes'])


    pos = nx.spring_layout(G, weight='unique_genes', seed=42)  # positions for all nodes

    # nodes
    nx.draw_networkx_nodes(G, pos, node_size=101, node_color='grey')

    nx.draw_networkx_edges(G, pos, edgelist=amr_only, width=4, edge_color='green')
    nx.draw_networkx_edges(G, pos, edgelist=vf_only, width=3, edge_color='orange')
    nx.draw_networkx_edges(G, pos, edgelist=both, width=3, edge_color="purple")
    nx.draw_networkx_edges(G, pos, edgelist=neither, width=2, edge_color='black', alpha=0.9)

    # labels
    #nx.draw_networkx_labels(G, pos, font_size=5, font_family="sans-serif")

    green_patch = mpatches.Patch(color='green', label='AMR in shared genes')
    orange_patch = mpatches.Patch(color='orange', label='VF in shared genes')
    purple_patch = mpatches.Patch(color='purple', label='AMR & VF in shared genes')
    black_patch = mpatches.Patch(color='black', label='Shared genes')

    plt.legend(handles=[green_patch, orange_patch, purple_patch, black_patch])

    plt.title('Orthologue Network')
    plt.xlabel('Edges weighted by number of unique shared orthologues')
    #plt.axis("off")

    nx.write_gpickle(G, 'network.pkl')
    nx.write_edgelist(G, 'network.edgelist')

    plt.savefig("network.png", dpi=300)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate highlighted orthology graph')
    parser.add_argument('--proteinortho_output', help='TSV output from proteinortho',
                        required=True)
    parser.add_argument('--rgi_output', help='TSV output from RGI',
                        required=True)
    parser.add_argument('--vf_blast_output', help='TSV output BLASTP vs VFDB',
                        required=True)
    parser.add_argument('--seq_metadata', default=None, help="TSV containing metadata (in '--metadata $column') to annotate each fasta ('name' column)")
    parser.add_argument('--metadata_col', default=None, help="Metadata column to use to annotate")


    args = parser.parse_args()

    if args.seq_metadata:
        if not args.metadata_col:
            print("\nIf --seq_metadata supplied then --metadata_col must be specified\n", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit(1)
    run(args)
