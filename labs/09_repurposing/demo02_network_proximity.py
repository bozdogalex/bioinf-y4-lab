"""
Demo 2 — Drug similarity network using Jaccard similarity on gene sets
"""

import pandas as pd
import networkx as nx
import itertools
import matplotlib.pyplot as plt

# Toy drug-gene mapping
drug2genes = {
    "Aspirin": {"PTGS1", "PTGS2", "IL6"},
    "Ibuprofen": {"PTGS2", "TNF"},
    "Statin": {"NPC1"}
}

# Build similarity network
G = nx.Graph()
G.add_nodes_from(drug2genes.keys())

for d1, d2 in itertools.combinations(drug2genes.keys(), 2):
    s1, s2 = drug2genes[d1], drug2genes[d2]
    jaccard = len(s1 & s2) / len(s1 | s2)
    if jaccard > 0:
        G.add_edge(d1, d2, weight=jaccard)

pos = nx.spring_layout(G, seed=42)
weights = [G[u][v]["weight"] * 5 for u, v in G.edges()]

nx.draw(G, pos, with_labels=True, width=weights, node_color="skyblue", node_size=800)
plt.title("Drug Similarity Network (Jaccard)")
plt.show()
