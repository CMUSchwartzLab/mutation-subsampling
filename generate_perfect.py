import os
import math
import random
import numpy as np
import argparse
from typing import List, Tuple, Dict

def build_full_binary_tree_canonical(n, seed):
    assert n >= 1
    
    rng = random.Random(seed)

    current = list(range(1, n+1))         # leaves
    next_internal = n + 1
    last_node = 2*n - 1

    parents = [0] * (last_node + 1)       # 1-based; 0 => no parent (root)
    edges: List[Tuple[int,int]] = []

    while len(current) > 1:
        a, b = rng.sample(current, 2)
        current.remove(a)
        current.remove(b)

        p = next_internal
        next_internal += 1
        edges.append((p, a))
        edges.append((p, b))
        parents[a] = p
        parents[b] = p
        current.append(p)

    root = current[0]
    if root != last_node:
        internals = sorted(set(v for e in edges for v in e) - set(range(1, n+1)))
        internals.remove(root)
        mapping = {root: last_node}
        target_ids = list(range(n+1, last_node))
        for old, new in zip(internals, target_ids):
            mapping[old] = new
        for leaf in range(1, n+1):
            mapping[leaf] = leaf

        new_edges = []
        new_parents = [0]*(last_node+1)
        for (u, v) in edges:
            u2, v2 = mapping[u], mapping[v]
            new_edges.append((u2, v2))
            new_parents[v2] = u2
        edges = new_edges
        parents = new_parents
        parents[last_node] = 0

    leaves = list(range(1, n+1))
    return parents, edges, leaves

def adjacency_matrix_E(n, edges):
    N = 2*n - 1
    E = np.zeros((N, N), dtype=int)
    for (u, v) in edges:
        E[u-1, v-1] = 1
    return E

def ancestor_matrix_A(n, parents):
    N = 2*n - 1
    A = np.zeros((N, N), dtype=int)
    for v in range(1, N+1):
        u = parents[v]
        while u != 0:
            A[u-1, v-1] = 1
            u = parents[u]
    return A

def assign_mutations_once(n, edges, m, seed):
    rng = random.Random(seed)
    non_root_children = [v for (_, v) in edges]
    edge_mutations: Dict[int, List[int]] = {v: [] for v in non_root_children}
    
    for j in range(m):
        v = rng.choice(non_root_children)
        edge_mutations[v].append(j)
    return edge_mutations

def build_mutation_matrix_all_nodes(n, edges, edge_mutations, m):
    N = 2*n - 1
    X = np.zeros((N, m), dtype=int)
    children = [[] for _ in range(N+1)]
    for (u, v) in edges:
        children[u].append(v)
    for v, muts in edge_mutations.items():
        if not muts:
            continue
        stack = [v]
        seen = []
        while stack:
            x = stack.pop()
            seen.append(x)
            stack.extend(children[x])
        rows = [idx-1 for idx in seen]
        X[np.ix_(rows, muts)] = 1
    return X

def leaf_only_matrix(n, X_all):
    return X_all[:n, :]

def export_dot_numeric(edges, edge_mutations, n):
    N = 2*n - 1
    lines = ["digraph G {"]
    for v in range(1, N+1):
        lines.append(f'  "{v}";')
    for (u, v) in edges:
        lab = ",".join(str(m) for m in edge_mutations.get(v, []))
        if lab:
            lines.append(f'  "{u}" -> "{v}" [label="{lab}"];')
        else:
            lines.append(f'  "{u}" -> "{v}";')
    lines.append("}")
    return "\n".join(lines)



def harmonic_number(k):
    return sum(1.0/i for i in range(1, k+1)) if k > 0 else 0.0

def run_batch():
    LEAF_LIST = [2, 3, 4, 5, 6, 10, 15, 20, 25, 30, 40, 50]
    OUT_ROOT = "../simulation_data/CLONE/"
    N_SIMS = 10

    BASE_SEED_TREE = 1337
    BASE_SEED_MUTS = 4242
    for leaf in LEAF_LIST:
        total_nodes = 2*leaf - 1
        n_edge = total_nodes - 1                 # = 2*leaf - 2
        m = int(math.ceil(n_edge * harmonic_number(n_edge)))
        clone_dir = os.path.join(OUT_ROOT, f"clone_{total_nodes}")
        os.makedirs(clone_dir, exist_ok=True)

        for sim in range(1, N_SIMS + 1):
            sim_dir = os.path.join(clone_dir, f"sim{sim}")
            os.makedirs(sim_dir, exist_ok=True)

            seed_tree = BASE_SEED_TREE + 1000*leaf + sim
            seed_muts = BASE_SEED_MUTS + 2000*leaf + sim
            
            parents, edges, leaves = build_full_binary_tree_canonical(leaf, seed=seed_tree) # Build tree

            edge_muts = assign_mutations_once(leaf, edges, m, seed=seed_muts) # Assign mutations exactly once

            X_all = build_mutation_matrix_all_nodes(leaf, edges, edge_muts, m) # Build all-nodes perfect-phylogeny matrix

            dot_str = export_dot_numeric(edges, edge_muts, leaf) # 4) Save DOT
            dot_path = os.path.join(sim_dir, "perfect_tree.dot")
            with open(dot_path, "w") as f:
                f.write(dot_str)

            mat_path = os.path.join(sim_dir, "perfect_matrix.tsv")
            np.savetxt(mat_path, X_all, fmt="%d", delimiter="\t")

            print(f"[leaf={leaf} | nodes={total_nodes} | sim={sim}]  m={m}  ->  {dot_path} , {mat_path}")

def parse_args():
    parser = argparse.ArgumentParser(description="Generate a perfect phylogeny and mutation matrix.")

    parser.add_argument("-n", type=int, required=True,
                        help="Number of leaves in the binary tree.")
    parser.add_argument("-m", type=int, required=True,
                        help="Number of mutations to assign.")
    parser.add_argument("-t", "--seed_tree", type=int, default=7,
                        help="Random seed for tree generation.")
    parser.add_argument("-s", "--seed_muts", type=int, default=11,
                        help="Random seed for mutation assignment.")
    parser.add_argument("-o", "--out_prefix", type=str, default="perfect",
                        help="Prefix for output files (DOT + TSV).")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    n = args.n
    m = args.m
    seed_tree = args.seed_tree
    seed_muts = args.seed_muts
    out_prefix = args.out_prefix

    parents, edges, leaves = build_full_binary_tree_canonical(n, seed_tree)

    edge_muts = assign_mutations_once(n, edges, m, seed_muts)

    X_all = build_mutation_matrix_all_nodes(n, edges, edge_muts, m)

    print("Edges:", edges)
    print("Mutation matrix shape:", X_all.shape)
    
    dot_str = export_dot_numeric(edges, edge_muts, n)
    dot_path = f"{out_prefix}_tree.dot"
    with open(dot_path, "w") as f:
        f.write(dot_str)
    
    mat_path = f"{out_prefix}_matrix.tsv"
    np.savetxt(mat_path, X_all, fmt="%d", delimiter="\t")

    print(f"Saved DOT to: {dot_path}")
    print(f"Saved matrix to: {mat_path}")
    