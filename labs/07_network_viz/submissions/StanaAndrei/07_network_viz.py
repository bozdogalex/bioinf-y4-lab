from __future__ import annotations
import gc
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# --------------------------
# CONFIGURATION
# --------------------------
HANDLE = "StanaAndrei"

# INPUTS (Must match where Lab 6 saved them)
# Note: Adjusted paths to match your previous Lab 6 structure
EXPR_CSV = Path(f"data/work/{HANDLE}/lab06/gene_expr_matrix.csv")
MODULES_CSV = Path(f"labs/06_networks/submissions/{HANDLE}/modules_{HANDLE}.csv")

# OUTPUTS
OUT_DIR = Path(f"labs/07_networkviz/submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"

# VISUALIZATION SETTINGS
CORR_METHOD = "spearman"
ADJ_THRESHOLD = 0.75       # Keep high to see clear structure
TOPK_HUBS = 10             # How many top genes to label
MAX_NODES_TO_PLOT = 1500   # Limit for plotting (prevents "hairball" and freezing)
SEED = 42

# --------------------------
# HELPER FUNCTIONS
# --------------------------

def read_modules(path: Path) -> Dict[str, int]:
    """Reads the gene-to-module mapping first."""
    if not path.exists():
        raise FileNotFoundError(f"Modules file not found: {path}\nRun Lab 6 first!")
    
    df = pd.read_csv(path)
    # Ensure columns exist
    if "Gene" not in df.columns or "Module" not in df.columns:
        # Fallback if headers are different, assume col 0 is Gene, col 1 is Module
        return dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(int)))
        
    return dict(zip(df["Gene"].astype(str), df["Module"].astype(int)))

def load_filtered_expression(path: Path, target_genes: Set[str]) -> pd.DataFrame:
    """
    Reads the huge CSV but ONLY keeps genes that are in our modules.
    This prevents Memory Errors.
    """
    print(f"   Reading expression data filtered by {len(target_genes)} module genes...")
    chunks = []
    # Read in chunks
    for chunk in pd.read_csv(path, index_col=0, chunksize=2000):
        # Find intersection
        common = chunk.index.intersection(target_genes)
        if not common.empty:
            subset = chunk.loc[common].astype(np.float32)
            # Log transform immediately
            np.log2(subset + 1, out=subset)
            chunks.append(subset)
    
    if not chunks:
        raise ValueError("No matching genes found in expression matrix!")
        
    df = pd.concat(chunks)
    print(f"   Loaded Matrix Shape: {df.shape}")
    return df

def build_graph(expr: pd.DataFrame, threshold: float) -> nx.Graph:
    """
    Calculates correlation and builds graph.
    Since 'expr' is already filtered to rank ~1000-3000 genes, 
    standard matrix calc is safe.
    """
    print(f"   Calculating {CORR_METHOD} correlation...")
    corr = expr.T.corr(method=CORR_METHOD)
    
    # Use abs value
    matrix = corr.abs().values
    
    # Mask lower triangle to avoid duplicates
    mask = np.triu(np.ones(matrix.shape, dtype=bool), k=1)
    
    # Find edges
    rows, cols = np.where((matrix >= threshold) & mask)
    
    print(f"   Found {len(rows)} strong connections.")
    
    G = nx.Graph()
    genes = corr.index.to_numpy()
    G.add_nodes_from(genes)
    G.add_edges_from(zip(genes[rows], genes[cols]))
    
    # Remove isolated nodes (they clutter the visualization)
    G.remove_nodes_from(list(nx.isolates(G)))
    
    return G

def get_hubs(G: nx.Graph, top_k: int) -> pd.DataFrame:
    """Identifies top connected genes."""
    degrees = dict(G.degree())
    # Sort by degree
    sorted_genes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
    
    # Create DataFrame
    df = pd.DataFrame(sorted_genes, columns=["Gene", "Degree"])
    return df.head(top_k)

# --------------------------
# MAIN EXECUTION
# --------------------------
if __name__ == "__main__":
    gc.enable()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    try:
        # 1. Load Modules
        print("1. Loading Module Definitions...")
        gene2module = read_modules(MODULES_CSV)
        target_genes = set(gene2module.keys())
        print(f"   Found {len(target_genes)} genes assigned to modules.")

        # 2. Load Expression Data (Filtered)
        print("2. Loading Expression Data...")
        expr_data = load_filtered_expression(EXPR_CSV, target_genes)
        
        # 3. Build Graph
        print("3. Building Network...")
        G = build_graph(expr_data, ADJ_THRESHOLD)
        print(f"   Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")

        # FREE MEMORY
        del expr_data
        gc.collect()

        # 4. Visualization Prep
        if G.number_of_nodes() == 0:
            print("   ERROR: Graph is empty. Try lowering ADJ_THRESHOLD.")
            sys.exit(0)

        # OPTIMIZATION: If graph is huge, keep only the strongest nodes for the plot
        if G.number_of_nodes() > MAX_NODES_TO_PLOT:
            print(f"   Graph too large for clean plot ({G.number_of_nodes()} nodes).")
            print(f"   Subsampling top {MAX_NODES_TO_PLOT} connected nodes...")
            degrees = dict(G.degree())
            top_nodes = sorted(degrees, key=degrees.get, reverse=True)[:MAX_NODES_TO_PLOT]
            G = G.subgraph(top_nodes).copy()

        # 5. Determine Colors and Sizes
        # Color map: Tab10 or Set3
        cmap = plt.get_cmap("tab10")
        
        node_colors = []
        node_sizes = []
        hubs_df = get_hubs(G, TOPK_HUBS)
        hub_set = set(hubs_df["Gene"])
        
        for node in G.nodes():
            mod_id = gene2module.get(node, 0)
            # Color based on module ID (modulo 10 to fit colormap)
            color = cmap((mod_id - 1) % 10)
            node_colors.append(color)
            
            # Size: Bigger if Hub
            if node in hub_set:
                node_sizes.append(150) # Hub size
            else:
                node_sizes.append(30)  # Standard size

        # 6. Plotting
        print("4. Calculating Layout (this may take a moment)...")
        plt.figure(figsize=(14, 12))
        
        # Spring layout is standard for organic looking networks
        pos = nx.spring_layout(G, seed=SEED, k=0.15, iterations=50)
        
        # Draw edges (thin and transparent)
        nx.draw_networkx_edges(G, pos, alpha=0.1, width=0.5, edge_color="gray")
        
        # Draw nodes
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, linewidths=0.5, edgecolors="white")
        
        # Label ONLY hubs (otherwise text is unreadable)
        hub_labels = {node: node for node in G.nodes() if node in hub_set}
        
        # Use a background box for labels so they are readable
        text_items = nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=9, font_weight="bold")
        for t in text_items.values():
            t.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

        plt.title(f"Gene Co-expression Network ({len(set(gene2module.values()))} Modules)\nTop {TOPK_HUBS} Hubs Labeled", fontsize=14)
        plt.axis("off")
        
        plt.tight_layout()
        plt.savefig(OUT_PNG, dpi=300)
        plt.close()
        print(f"SUCCESS! Network image saved to: {OUT_PNG}")
        
        # 7. Save Hubs
        hubs_df.to_csv(OUT_HUBS, index=False)
        print(f"Hub list saved to: {OUT_HUBS}")

    except Exception as e:
        print(f"CRITICAL ERROR: {e}")
        import traceback
        traceback.print_exc()