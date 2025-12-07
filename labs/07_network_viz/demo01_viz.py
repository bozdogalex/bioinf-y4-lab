"""
Demo 1 — Vizualizarea unei rețele toy cu module și hub genes

Scop:
- Arată cum să colorăm nodurile după modul
- Cum să identificăm genele hub (grad mare)
- Cum să exportăm figura

Acest demo folosește un graf mic, nu date reale.
"""

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

# 1) Construim un graf mic
G = nx.Graph()
edges = [
    ("GeneA", "GeneB"),
    ("GeneA", "GeneC"),
    ("GeneB", "GeneC"),
    ("GeneC", "GeneD"),
    ("GeneD", "GeneE"),
    ("GeneE", "GeneF"),
]
G.add_edges_from(edges)

# 2) Module asignate manual
gene2module = {
    "GeneA": 1, "GeneB": 1, "GeneC": 1,
    "GeneD": 2, "GeneE": 2, "GeneF": 2,
}

# 3) Culori după modul
cmap = plt.get_cmap("tab10")
node_colors = [cmap((gene2module[n] - 1) % 10) for n in G.nodes()]

# 4) Hub genes (grad mare)
deg = dict(G.degree())
hubs = sorted(deg.items(), key=lambda x: x[1], reverse=True)[:2]  # top 2
hub_nodes = {n for n, _ in hubs}
print("Hub genes:", hubs)

# 5) Layout + vizualizare
pos = nx.spring_layout(G, seed=42)
plt.figure(figsize=(6, 5))
nx.draw_networkx_edges(G, pos, alpha=0.3)
nx.draw_networkx_nodes(
    G, pos,
    node_color=node_colors,
    node_size=[300 if n in hub_nodes else 150 for n in G.nodes()]
)
nx.draw_networkx_labels(G, pos, font_size=10)

plt.axis("off")
plt.tight_layout()
plt.savefig("demo_network.png", dpi=200)
plt.show()

# 6) Export hub genes în CSV
df_hubs = pd.DataFrame(hubs, columns=["Gene", "Degree"])
df_hubs.to_csv("demo_hubs.csv", index=False)
print("Am salvat demo_network.png și demo_hubs.csv")
