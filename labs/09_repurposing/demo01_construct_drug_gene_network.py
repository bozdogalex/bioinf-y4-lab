"""
Demo 1 — Construirea unei rețele bipartite drug–gene (toy dataset)
"""

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Toy dataset
data = {
    "drug": ["Aspirin", "Aspirin", "Aspirin", "Ibuprofen", "Ibuprofen", "Statin"],
    "gene": ["PTGS1", "PTGS2", "IL6", "PTGS2", "TNF", "NPC1"]
}
df = pd.DataFrame(data)

# Build bipartite graph
G = nx.Graph()
G.add_nodes_from(df["drug"], bipartite="drug")
G.add_nodes_from(df["gene"], bipartite="gene")
G.add_edges_from(zip(df["drug"], df["gene"]))

# Draw
plt.figure(figsize=(7, 5))
pos = nx.spring_layout(G, seed=42)

nx.draw(
    G, pos,
    node_color=["skyblue" if n in df["drug"].unique() else "salmon" for n in G.nodes()],
    node_size=600, with_labels=True
)
plt.title("Bipartite Drug–Gene Network (Toy Example)")
plt.show()
