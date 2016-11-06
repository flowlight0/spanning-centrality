import networkx as nx
import scipy
import numpy as np

lines = open('karate.txt').readlines()
es = [list(map(int, line.strip().split('\t'))) for line in lines]

def num_spanning_trees(vs, es):
    g = nx.Graph()
    g.add_nodes_from(vs)
    for u, v in es:
        g.add_edge(u, v)
    l = nx.laplacian_matrix(g)[1:, 1:].toarray()
    return np.linalg.det(l)

vs = set()
for a, b in es:
    vs.add(a)
    vs.add(b)

all = num_spanning_trees(vs, es)

for a, b in es:
    es_ = [[u, v] for u, v in es if u != a or v != b]
    print("{0}\t{1}\t{2:.3f}".format(a, b, 1 - num_spanning_trees(vs, es_) / all))
