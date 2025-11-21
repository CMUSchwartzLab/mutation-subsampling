from __future__ import annotations
import argparse
import random
from typing import Dict, List, Tuple, Optional, Set
import pydot


def read_matrix_tsv(path):
    D: List[List[int]] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            row = [int(x) for x in line.split()]
            D.append(row)
    if not D:
        raise ValueError("Empty matrix file.")
    
    m = len(D[0])
    for r, row in enumerate(D):
        if len(row) != m:
            raise ValueError(f"Non-rectangular matrix at row {r}: expected {m} cols, got {len(row)}.")
    return D  # R x M



def write_matrix_tsv(path, D):
    with open(path, "w") as f:
        for row in D:
            f.write("\t".join(str(x) for x in row) + "\n")


# ---------------- DOT helpers ----------------

def _clean(s):
    return "" if s is None else s.strip().strip('"')


def _parse_edge_label_to_ints(lbl):
    if not lbl:
        return []
    s = _clean(lbl)
    if not s:
        return []
    toks = s.replace("\\n", "\n").replace(",", " ").split()
    out: List[int] = []
    for tok in toks:
        num = ""
        for ch in tok:
            if ch.isdigit() or (ch == "-" and not num):
                num += ch
            else:
                break
        if num:
            out.append(int(num))
    return out


def _format_edge_label_with_base(gains0, losses0, base):
    parts= []
    for c0 in sorted(gains0):
        parts.append(f"{c0 + base}g")
    for c0, lidx in sorted(losses0):
        parts.append(f"{c0 + base}l{lidx}")
    return "\\n".join(parts)


class DotTree:
    def __init__(self):
        self.children: Dict[int, List[int]] = {}
        self.parent: Dict[int, int] = {}
        self.root: int = -1
        self.edge_gains_raw: Dict[Tuple[int, int], List[int]] = {}

    @staticmethod
    def from_dot_file(dot_path: str) -> "DotTree":
        graphs = pydot.graph_from_dot_file(dot_path)
        if not graphs:
            raise RuntimeError(f"Could not read DOT file: {dot_path}")
        g: pydot.Dot = graphs[0]
        t = DotTree()

        node_ids: Set[int] = set()
        for node in g.get_nodes():
            nm = _clean(node.get_name())
            if nm in ("graph", "node", "edge"):
                continue
            try:
                vid = int(nm)
            except ValueError:
                raise ValueError(f"DOT node name is not an integer id: {nm}")
            node_ids.add(vid)

        indeg: Dict[int, int] = {vid: 0 for vid in node_ids}

        for edge in g.get_edges():
            u = int(_clean(edge.get_source()))
            v = int(_clean(edge.get_destination()))
            t.children.setdefault(u, []).append(v)
            t.parent[v] = u
            indeg[v] = indeg.get(v, 0) + 1
            indeg.setdefault(u, 0)
            t.edge_gains_raw[(u, v)] = _parse_edge_label_to_ints(edge.get_label())
            node_ids.add(u)
            node_ids.add(v)

        for vid in list(node_ids):
            t.children.setdefault(vid, [])
            indeg.setdefault(vid, 0)

        roots = [vid for vid, d in indeg.items() if d == 0]
        if len(roots) != 1:
            raise ValueError(f"DOT must have exactly one root (found roots={roots})")
        t.root = roots[0]
        return t

    def leaves(self):
        return [v for v, ch in self.children.items() if len(ch) == 0]



class DolloStrictFull:
    def __init__(self, D_perfect, tree, mut_base):
        self.D = D_perfect              # R x M (0/1)
        self.R = len(D_perfect)         # nodes count
        self.M = len(D_perfect[0])      # mutations count
        self.T = tree
        if mut_base not in (0, 1):
            raise ValueError("--mut-base must be 0 or 1")
        self.mut_base = mut_base

        # STRICT node set: must be exactly {1..R}
        node_ids = set(self.T.children.keys())
        #print(node_ids)
        for kids in self.T.children.values():
            node_ids.update(kids)
        expected_nodes = set(range(1, self.R + 1))
        if node_ids != expected_nodes:
            raise ValueError(f"DOT node IDs must be exactly 1..R. Found={sorted(node_ids)}, expected={sorted(expected_nodes)}")

        # Leaves must be 1..L where R = 2L - 1
        if (self.R + 1) % 2 != 0:
            raise ValueError(f"Row count R={self.R} is not of the form 2L-1.")
        self.L = (self.R + 1) // 2
        leaves = sorted(self.T.leaves())
        expected_leaves = list(range(1, self.L + 1))
        if leaves != expected_leaves:
            raise ValueError(f"Leaves must be 1..L with L={self.L}. Found leaves={leaves}")

        if self.T.root != self.R:
            raise ValueError(f"Root must be node {self.R} (last row). Found root={self.T.root}")

        
        self.edge_gains: Dict[Tuple[int, int], List[int]] = {}
        for e, raw in self.T.edge_gains_raw.items():
            norm: List[int] = []
            for C in raw:
                c0 = C - self.mut_base
                if not (0 <= c0 < self.M):
                    raise ValueError(
                        f"Gain index {C} on edge {e} invalid for mut_base={self.mut_base} "
                        f"(normalized {c0} not in [0..{self.M-1}])."
                    )
                norm.append(c0)
            self.edge_gains[e] = norm
        self.a: Dict[int, List[int]] = {}  # {node -> [M ints in {0,1,2..k+1}]}
        self.b: Dict[int, List[int]] = {}  # {node -> [M ints in {0,1}]}
        self.edge_losses: Dict[Tuple[int, int], List[Tuple[int, int]]] = {}  # (c0, loss_idx)

    def _init_states_from_gains(self):
        self.a[self.T.root] = [0] * self.M
        self.b[self.T.root] = [0] * self.M
        stack = [self.T.root]
        while stack:
            u = stack.pop()
            for v in self.T.children[u]:
                self.a[v] = list(self.a[u])
                self.b[v] = list(self.b[u])
                for c0 in self.edge_gains.get((u, v), []):
                    self.a[v][c0] = 1
                    self.b[v][c0] = 1
                stack.append(v)

    def simulate_losses(self, k: int, loss_rate: float, seed: int) -> None:
        rng = random.Random(seed)
        remaining = [k] * self.M
        used = [0] * self.M

        def dfs(u: int):
            for v in self.T.children[u]:
                #print(v)
                # no losses on leaves (v is leaf if it has no children)
                if len(self.T.children[v]) == 0:
                    continue

                gains_uv = set(self.edge_gains.get((u, v), []))
                #print(u, v, gains_uv)
                present = [c0 for c0 in range(self.M)
                           if self.a[v][c0] == 1 and c0 not in gains_uv and remaining[c0] > 0]
                #print("=> ", present, self.M, self.a[v][:],gains_uv )
                for c0 in present:
                    #print("here")
                    if rng.random() <= loss_rate:
                        used[c0] += 1
                        if used[c0] <= k:
                            self.edge_losses.setdefault((u, v), []).append((c0, used[c0]))
                            remaining[c0] -= 1

                # apply losses to child v then push down
                for (c0, lidx) in self.edge_losses.get((u, v), []):
                    self.a[v][c0] = 1 + lidx    # 2..k+1
                    self.b[v][c0] = 0

                self._push_down(v)
                dfs(v)

        self._init_states_from_gains()
        dfs(self.T.root)

    def _push_down(self, start):
        stack = [start]
        while stack:
            u = stack.pop()
            for v in self.T.children[u]:
                self.a[v] = list(self.a[u])
                self.b[v] = list(self.b[u])
                for c0 in self.edge_gains.get((u, v), []):
                    self.a[v][c0] = 1
                    self.b[v][c0] = 1
                for (c0, lidx) in self.edge_losses.get((u, v), []):
                    self.a[v][c0] = 1 + lidx
                    self.b[v][c0] = 0
                stack.append(v)

    
    def make_A_full(self):
        A = [[0] * self.M for _ in range(self.R)]
        for v in range(1, self.R + 1):
            if v in self.a:
                A[v - 1] = list(self.a[v])
        return A

    def make_B_full(self):
        B = [[0] * self.M for _ in range(self.R)]
        for v in range(1, self.R + 1):
            if v in self.b:
                B[v - 1] = list(self.b[v])
        return B

    def write_dollo_dot(self, out_dot):
        graph = pydot.Dot(graph_type="digraph", strict=True)
        for v in range(1, self.R + 1):
            graph.add_node(pydot.Node(str(v)))
        for u, kids in self.T.children.items():
            for v in kids:
                gains0 = self.edge_gains.get((u, v), [])
                losses0 = self.edge_losses.get((u, v), [])
                label = _format_edge_label_with_base(gains0, losses0, self.mut_base)
                if label:
                    graph.add_edge(pydot.Edge(str(u), str(v), label=label))
                else:
                    graph.add_edge(pydot.Edge(str(u), str(v)))
        graph.write_raw(out_dot)


def main():
    ap = argparse.ArgumentParser(description="STRICT k-Dollo simulation (full-node matrices).")
    ap.add_argument("perfect_matrix_tsv", help="TSV, shape (2L-1) x M, 0/1 entries, no header/index")
    ap.add_argument("perfect_dot", help="DOT with node IDs 1..(2L-1), edges labeled with integer mutation indices")
    ap.add_argument("-k", "--k", type=int, required=True, help="Max losses per character")
    ap.add_argument("--loss", type=float, required=True, help="Loss rate per eligible internal edge (0..1)")
    ap.add_argument("-s", "--seed", type=int, default=0, help="Random seed")
    ap.add_argument("--mut-base", type=int, choices=[0, 1], default=1,
                    help="Base for mutation labels in DOT (0 or 1). Default: 1")
    ap.add_argument("-A", "--outA", required=True, help="Output A_full TSV (2L-1 x M, entries 0/1/2..k+1)")
    ap.add_argument("-B", "--outB", default="", help="Optional output B_full TSV (2L-1 x M, binary)")
    ap.add_argument("--dot", required=True, help="Output Dollo DOT (labels Cg / Cl#; C printed with --mut-base)")
    args = ap.parse_args()

    D = read_matrix_tsv(args.perfect_matrix_tsv)  # R x M
    T = DotTree.from_dot_file(args.perfect_dot)

    sim = DolloStrictFull(D, T, mut_base=args.mut_base)
    sim.simulate_losses(k=args.k, loss_rate=args.loss, seed=args.seed)

    A_full = sim.make_A_full()
    write_matrix_tsv(args.outA, A_full)

    if args.outB:
        B_full = sim.make_B_full()
        write_matrix_tsv(args.outB, B_full)

    sim.write_dollo_dot(args.dot)

    print(f"Wrote A_full: {args.outA}")
    if args.outB:
        print(f"Wrote B_full: {args.outB}")
    print(f"Wrote DOT: {args.dot}")


if __name__ == "__main__":
    main()
