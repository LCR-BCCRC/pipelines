from phyclone.map import get_map_node_ccfs
from phyclone.run import get_labels_table
from phyclone.tree import Tree

import gzip
import networkx as nx
import numpy as np
import pandas as pd
import pickle


def main(args):
    with gzip.GzipFile(args.in_file, "rb") as fh:
        results = pickle.load(fh)

    out_df = []

    for state in results["trace"][args.burnin:]:
        phyclone_tree = Tree.from_dict(results["data"], state["tree"])

        tree = get_tree_from_phyclone_tree(phyclone_tree)

        prevs = np.array([tree.nodes[n]["clonal_prev"][0] for n in tree.nodes])

        prevs += 1e-6

        prevs = prevs / np.sum(prevs)

        entropy = -np.sum(prevs * np.log2(prevs))

        num_clones = len(prevs)

        labels = get_labels_table(results["data"], phyclone_tree, clusters=results["clusters"])

        num_snvs = labels.groupby("cluster_id")["mutation_id"].nunique().values

        min_num_snvs = min(num_snvs)

        max_num_snvs = max(num_snvs)

        mean_num_snvs = np.mean(num_snvs)

        median_num_snvs = np.median(num_snvs)

        out_df.append({
            "entropy": entropy,
            "num_clones": num_clones,
            "min_num_snvs": min_num_snvs,
            "max_num_snvs": max_num_snvs,
            "mean_num_snvs": mean_num_snvs,
            "median_num_snvs": median_num_snvs
        })

    out_df = pd.DataFrame(out_df)

    out_df = pd.DataFrame([out_df.mean()])

    out_df.insert(0, "patient_id", args.patient_id)

    out_df.to_csv(args.out_file, index=False, sep="\t")


def get_tree_from_phyclone_tree(tree):
    """ Convert a Phyclone tree object to a graph for benchmarking.

    Parameters
    ----------
    tree: (phyclone.tree.Tree)

    Returns
    -------
    nx.Digraph representing clone phylogeny with nodes "snvs", "cellular_prev" and "clonal_prev" set for each node
    """
    ccfs = get_map_node_ccfs(tree)

    G = nx.DiGraph()

    for n in tree.graph.nodes:
        G.add_node(n)

        G.nodes[n]["cellular_prev"] = ccfs[n]

        G.nodes[n]["snvs"] = [x.name for x in tree.node_data[n]]

    for u, v in tree.graph.edges:
        G.add_edge(u, v)

    roots = []

    for n in G.nodes:
        if G.in_degree(n) == 0:
            roots.append(n)

    for r in roots:
        set_clonal_prev(G, r)

    return G


def set_clonal_prev(G, node):
    clonal_prev = G.nodes[node]["cellular_prev"].copy()

    for child in G.successors(node):
        clonal_prev -= G.nodes[child]["cellular_prev"]

        set_clonal_prev(G, child)

    G.nodes[node]["clonal_prev"] = clonal_prev


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-p", "--patient-id", required=True)

    parser.add_argument("-b", "--burnin", default=0, type=int)

    cli_args = parser.parse_args()

    main(cli_args)
