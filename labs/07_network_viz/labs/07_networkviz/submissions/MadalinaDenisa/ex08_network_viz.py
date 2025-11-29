from pathlib import Path
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

HANDLE = "MadalinaDenisa"

# Input files
EXPR_CSV = Path("expression_matrix.csv")
MODULES_CSV = Path("modules_MadalinaDenisa.csv")

# Adjacency parameters
CORR_METHOD = "spearman"
USE_ABS_CORR = True
ADJ_THRESHOLD = 0.2  # prag mai permisiv
WEIGHTED = False

# Visualization
SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.15
OUT_PNG = Path(f"network_{HANDLE}.png")
OUT_HUBS = Path(f"hubs_{HANDLE}.csv")

def ensure_exists(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

def read_expression_matrix(path):
    df = pd.read_csv(path, index_col=0)
    return df

def read_modules_csv(path):
    df = pd.read_csv(path)
    return dict(zip(df['Gene'], df['Module']))

def correlation_to_adjacency(expr, threshold=0.2, use_abs=True, weighted=False, method="spearman"):
    if method == "spearman":
        corr = expr.T.corr(method='spearman')
    else:
        corr = expr.T.corr()
    if use_abs:
        corr = corr.abs()
    if not weighted:
        adj = (corr >= threshold).astype(int)
    else:
        adj = corr.where(corr >= threshold, 0)
    np.fill_diagonal(adj.values, 0)
    return adj

def graph_from_adjacency(adj):
    G = nx.from_pandas_adjacency(adj)
    G.remove_nodes_from(list(nx.isolates(G)))
    return G

def color_map_from_modules(nodes, gene2module):
    import matplotlib.cm as cm
    cmap = cm.get_cmap('tab10')
    colors = {}
    module_ids = list(set(gene2module.values()))
    module_color = {m: cmap(i/len(module_ids)) for i, m in enumerate(module_ids)}
    for n in nodes:
        colors[n] = module_color.get(gene2module.get(n, -1), (0.5,0.5,0.5))
    return colors

def compute_hubs(G, topk=10):
    deg = dict(G.degree())
    df = pd.DataFrame(deg.items(), columns=["Gene","Degree"]).sort_values("Degree", ascending=False)
    return df.head(topk)

# ---------------------
# Main
# ---------------------
ensure_exists(EXPR_CSV)
ensure_exists(MODULES_CSV)

expr = read_expression_matrix(EXPR_CSV)
gene2module = read_modules_csv(MODULES_CSV)
adj = correlation_to_adjacency(expr, threshold=ADJ_THRESHOLD, use_abs=USE_ABS_CORR, weighted=WEIGHTED, method=CORR_METHOD)
G = graph_from_adjacency(adj)

print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

node_colors = [color_map_from_modules(G.nodes(), gene2module)[n] for n in G.nodes()]
hubs_df = compute_hubs(G, TOPK_HUBS)
node_sizes = [NODE_BASE_SIZE*2 if n in hubs_df['Gene'].values else NODE_BASE_SIZE for n in G.nodes()]

plt.figure(figsize=(10,10))
pos = nx.spring_layout(G, seed=SEED)
nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA)
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes)
nx.draw_networkx_labels(G, pos, labels={n:n for n in hubs_df['Gene']}, font_size=10)
plt.axis('off')
plt.tight_layout()
plt.savefig(OUT_PNG)
hubs_df.to_csv(OUT_HUBS, index=False)
print(f"Network visualization exported: {OUT_PNG}")
print(f"Hub genes saved: {OUT_HUBS}")
