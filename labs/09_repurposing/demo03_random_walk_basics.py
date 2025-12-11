"""
Demo 3 — Simplified RWR for education
"""

import numpy as np
import networkx as nx

# Toy graph
G = nx.Graph()
G.add_edges_from([("Disease", "G1"), ("G1", "DrugA"), ("G1", "DrugB"), ("G2", "DrugC")])

nodes = list(G.nodes())
n = len(nodes)
idx = {n: i for i, n in enumerate(nodes)}

# Adjacency matrix
A = nx.to_numpy_array(G, nodelist=nodes)
A = A / A.sum(axis=1, keepdims=True)

# RWR parameters
alpha = 0.85
p = np.zeros(n)
p[idx["Disease"]] = 1.0

for _ in range(50):
    p = alpha * A.T.dot(p) + (1 - alpha) * (np.arange(n) == idx["Disease"])

ranking = sorted(zip(nodes, p), key=lambda x: x[1], reverse=True)
print("RWR Ranking:", ranking)
