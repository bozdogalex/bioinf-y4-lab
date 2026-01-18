"""
Bonus — Random Walk with Restart (RWR) on Drug-Gene Network

Implementează RWR pentru ranking-ul medicamentelor și compară cu metoda de proximitate.

Algorithm:
- Start from disease genes (seed nodes)
- Random walk with probability α of restarting to seed nodes
- Compute steady-state probability for each node
- Rank drugs by their RWR score
"""

from __future__ import annotations
from pathlib import Path
from typing import Set, Dict, Tuple
import pickle

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import gmres

# --------------------------
# Configuration
# --------------------------
HANDLE = "AlexTGoCreative"

# Input files
GRAPH_PICKLE = Path(f"labs/09_repurposing/submissions/{HANDLE}/assignment/network_drug_gene_{HANDLE}.gpickle")
DISEASE_GENES_TXT = Path(f"data/work/{HANDLE}/lab09/disease_genes_{HANDLE}.txt")
DRUG_PRIORITY_CSV = Path(f"labs/09_repurposing/submissions/{HANDLE}/assignment/drug_priority_{HANDLE}.csv")

# Output directory
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}/assignment")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Output files
OUT_RWR_RANKING = OUT_DIR / f"bonus_rwr_ranking_{HANDLE}.csv"
OUT_COMPARISON = OUT_DIR / f"bonus_comparison_{HANDLE}.csv"
OUT_COMPARISON_PLOT = OUT_DIR / f"bonus_comparison_{HANDLE}.png"

# RWR Parameters
RESTART_PROB = 0.7  # Probability of restart (α)
MAX_ITER = 100      # Maximum iterations
TOLERANCE = 1e-6    # Convergence tolerance


# ========================================
# Load Data
# ========================================

def load_graph(path: Path) -> nx.Graph:
    """Load the bipartite drug-gene graph."""
    with open(path, 'rb') as f:
        return pickle.load(f)


def load_disease_genes(path: Path) -> Set[str]:
    """Load disease gene list."""
    genes = set()
    with open(path, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                genes.add(gene)
    return genes


def get_drug_nodes(G: nx.Graph) -> list:
    """Extract drug nodes from bipartite graph."""
    return [n for n, d in G.nodes(data=True) if d.get("bipartite") == "drug"]


# ========================================
# Random Walk with Restart Implementation
# ========================================

def create_transition_matrix(G: nx.Graph) -> Tuple[csr_matrix, Dict[str, int], Dict[int, str]]:
    """
    Create column-stochastic transition matrix for random walk.
    
    Returns:
        - Transition matrix (sparse)
        - node_to_idx mapping
        - idx_to_node mapping
    """
    nodes = list(G.nodes())
    n = len(nodes)
    
    # Create mappings
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    idx_to_node = {i: node for i, node in enumerate(nodes)}
    
    # Build adjacency matrix
    row_ind = []
    col_ind = []
    data = []
    
    for node in nodes:
        neighbors = list(G.neighbors(node))
        degree = len(neighbors)
        
        if degree > 0:
            node_idx = node_to_idx[node]
            for neighbor in neighbors:
                neighbor_idx = node_to_idx[neighbor]
                # Column-stochastic: normalize by out-degree
                row_ind.append(neighbor_idx)
                col_ind.append(node_idx)
                data.append(1.0 / degree)
    
    # Create sparse transition matrix
    transition_matrix = csr_matrix((data, (row_ind, col_ind)), shape=(n, n))
    
    return transition_matrix, node_to_idx, idx_to_node


def create_restart_vector(
    disease_genes: Set[str],
    node_to_idx: Dict[str, int],
    G: nx.Graph
) -> np.ndarray:
    """
    Create restart probability vector.
    Uniform distribution over disease genes that exist in the graph.
    """
    n = len(node_to_idx)
    restart_vec = np.zeros(n)
    
    # Find disease genes in graph
    disease_genes_in_graph = [g for g in disease_genes if g in G.nodes()]
    
    if not disease_genes_in_graph:
        raise ValueError("No disease genes found in graph!")
    
    # Uniform distribution over disease genes
    for gene in disease_genes_in_graph:
        idx = node_to_idx[gene]
        restart_vec[idx] = 1.0 / len(disease_genes_in_graph)
    
    return restart_vec


def random_walk_with_restart(
    G: nx.Graph,
    disease_genes: Set[str],
    restart_prob: float = 0.7,
    max_iter: int = 100,
    tol: float = 1e-6
) -> Dict[str, float]:
    """
    Perform Random Walk with Restart (RWR).
    
    Algorithm:
        p^(t+1) = (1 - α) * W * p^(t) + α * r
        
    where:
        p^(t) = probability vector at iteration t
        W = column-stochastic transition matrix
        r = restart vector (uniform over disease genes)
        α = restart probability
    
    Returns:
        Dictionary mapping node -> steady-state probability
    """
    print(f"[RWR] Running with α={restart_prob}, max_iter={max_iter}, tol={tol}")
    
    # Create transition matrix
    W, node_to_idx, idx_to_node = create_transition_matrix(G)
    n = len(node_to_idx)
    
    # Create restart vector
    r = create_restart_vector(disease_genes, node_to_idx, G)
    
    # Initialize probability vector (start from restart distribution)
    p = r.copy()
    
    # Iterative method
    for iteration in range(max_iter):
        p_old = p.copy()
        
        # RWR update: p = (1-α)Wp + αr
        p = (1 - restart_prob) * W.dot(p) + restart_prob * r
        
        # Check convergence
        diff = np.linalg.norm(p - p_old, ord=1)
        
        if diff < tol:
            print(f"[RWR] Converged after {iteration + 1} iterations (diff={diff:.2e})")
            break
    else:
        print(f"[RWR] Max iterations reached (diff={diff:.2e})")
    
    # Convert to dictionary
    scores = {idx_to_node[i]: p[i] for i in range(n)}
    
    return scores


# ========================================
# Ranking and Comparison
# ========================================

def rank_drugs_by_rwr(
    G: nx.Graph,
    rwr_scores: Dict[str, float]
) -> pd.DataFrame:
    """Rank drugs by their RWR scores."""
    drug_nodes = get_drug_nodes(G)
    
    results = []
    for drug in drug_nodes:
        score = rwr_scores.get(drug, 0.0)
        results.append({'drug': drug, 'rwr_score': score})
    
    df = pd.DataFrame(results)
    df = df.sort_values('rwr_score', ascending=False)  # Higher score = better
    df['rwr_rank'] = range(1, len(df) + 1)
    
    return df


def compare_rankings(
    df_rwr: pd.DataFrame,
    df_proximity: pd.DataFrame
) -> pd.DataFrame:
    """Compare RWR ranking with proximity ranking."""
    
    # Merge dataframes
    df_comparison = df_rwr.merge(
        df_proximity[['drug', 'distance']],
        on='drug',
        how='inner'
    )
    
    # Add proximity rank (lower distance = better rank)
    df_comparison = df_comparison.sort_values('distance')
    df_comparison['proximity_rank'] = range(1, len(df_comparison) + 1)
    
    # Sort by RWR rank
    df_comparison = df_comparison.sort_values('rwr_rank')
    
    # Calculate rank difference
    df_comparison['rank_diff'] = abs(df_comparison['rwr_rank'] - df_comparison['proximity_rank'])
    
    # Add agreement metric
    df_comparison['agreement'] = df_comparison['rank_diff'] == 0
    
    return df_comparison


# ========================================
# Visualization
# ========================================

def visualize_comparison(df_comparison: pd.DataFrame, output_path: Path):
    """Create comprehensive comparison visualization."""
    
    fig = plt.figure(figsize=(16, 10))
    
    # 1. Scatter plot: RWR rank vs Proximity rank
    ax1 = plt.subplot(2, 2, 1)
    ax1.scatter(df_comparison['proximity_rank'], df_comparison['rwr_rank'], 
                alpha=0.6, s=100, c='steelblue', edgecolors='black')
    
    # Add diagonal line (perfect agreement)
    max_rank = max(df_comparison['proximity_rank'].max(), df_comparison['rwr_rank'].max())
    ax1.plot([1, max_rank], [1, max_rank], 'r--', alpha=0.5, label='Perfect agreement')
    
    ax1.set_xlabel('Proximity Rank', fontsize=12)
    ax1.set_ylabel('RWR Rank', fontsize=12)
    ax1.set_title('Ranking Comparison', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # 2. Top 10 comparison
    ax2 = plt.subplot(2, 2, 2)
    top_n = 10
    top_drugs = df_comparison.head(top_n)
    
    x = np.arange(len(top_drugs))
    width = 0.35
    
    ax2.barh(x - width/2, top_drugs['rwr_rank'], width, label='RWR Rank', color='#4ECDC4')
    ax2.barh(x + width/2, top_drugs['proximity_rank'], width, label='Proximity Rank', color='#FF6B6B')
    
    ax2.set_yticks(x)
    ax2.set_yticklabels(top_drugs['drug'], fontsize=9)
    ax2.invert_yaxis()
    ax2.set_xlabel('Rank Position', fontsize=12)
    ax2.set_title(f'Top {top_n} Drugs - Ranking Comparison', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(axis='x', alpha=0.3)
    
    # 3. Rank difference distribution
    ax3 = plt.subplot(2, 2, 3)
    ax3.hist(df_comparison['rank_diff'], bins=20, color='coral', alpha=0.7, edgecolor='black')
    ax3.set_xlabel('Rank Difference (|RWR - Proximity|)', fontsize=12)
    ax3.set_ylabel('Frequency', fontsize=12)
    ax3.set_title('Distribution of Rank Differences', fontsize=14, fontweight='bold')
    ax3.grid(axis='y', alpha=0.3)
    
    # Add statistics
    mean_diff = df_comparison['rank_diff'].mean()
    median_diff = df_comparison['rank_diff'].median()
    ax3.axvline(mean_diff, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_diff:.1f}')
    ax3.axvline(median_diff, color='blue', linestyle='--', linewidth=2, label=f'Median: {median_diff:.1f}')
    ax3.legend()
    
    # 4. Score vs Distance scatter
    ax4 = plt.subplot(2, 2, 4)
    
    # Normalize RWR scores for better visualization
    df_comparison['rwr_score_norm'] = df_comparison['rwr_score'] / df_comparison['rwr_score'].max()
    
    scatter = ax4.scatter(df_comparison['distance'], df_comparison['rwr_score_norm'],
                         c=df_comparison['rank_diff'], cmap='RdYlGn_r', 
                         s=100, alpha=0.7, edgecolors='black')
    
    ax4.set_xlabel('Proximity Distance (lower = better)', fontsize=12)
    ax4.set_ylabel('RWR Score (normalized, higher = better)', fontsize=12)
    ax4.set_title('Method Correlation', fontsize=14, fontweight='bold')
    ax4.grid(alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Rank Difference', fontsize=10)
    
    # Add annotations for top candidates
    top_5 = df_comparison.head(5)
    for _, row in top_5.iterrows():
        ax4.annotate(row['drug'], 
                    (row['distance'], row['rwr_score_norm']),
                    fontsize=8, alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"[INFO] Saved comparison visualization to {output_path}")
    plt.close()


# ========================================
# Analysis and Reporting
# ========================================

def analyze_comparison(df_comparison: pd.DataFrame):
    """Analyze and print comparison statistics."""
    
    print("\n" + "="*80)
    print("COMPARISON ANALYSIS: RWR vs Proximity")
    print("="*80)
    
    # Overall statistics
    print("\n[OVERALL STATISTICS]")
    print(f"Total drugs ranked: {len(df_comparison)}")
    print(f"Mean rank difference: {df_comparison['rank_diff'].mean():.2f}")
    print(f"Median rank difference: {df_comparison['rank_diff'].median():.1f}")
    print(f"Max rank difference: {df_comparison['rank_diff'].max():.0f}")
    print(f"Perfect agreements (same rank): {df_comparison['agreement'].sum()} ({df_comparison['agreement'].sum()/len(df_comparison)*100:.1f}%)")
    
    # Correlation
    from scipy.stats import spearmanr, kendalltau
    spearman_corr, spearman_p = spearmanr(df_comparison['rwr_rank'], df_comparison['proximity_rank'])
    kendall_corr, kendall_p = kendalltau(df_comparison['rwr_rank'], df_comparison['proximity_rank'])
    
    print(f"\n[RANK CORRELATION]")
    print(f"Spearman correlation: {spearman_corr:.3f} (p={spearman_p:.2e})")
    print(f"Kendall tau: {kendall_corr:.3f} (p={kendall_p:.2e})")
    
    # Top 10 overlap
    top_10_rwr = set(df_comparison.nsmallest(10, 'rwr_rank')['drug'])
    top_10_prox = set(df_comparison.nsmallest(10, 'proximity_rank')['drug'])
    overlap = len(top_10_rwr & top_10_prox)
    
    print(f"\n[TOP 10 OVERLAP]")
    print(f"Drugs in both top 10: {overlap}/10 ({overlap/10*100:.0f}%)")
    print(f"RWR only: {top_10_rwr - top_10_prox}")
    print(f"Proximity only: {top_10_prox - top_10_rwr}")
    
    # Top 5 from each method
    print(f"\n[TOP 5 CANDIDATES]")
    print("\nRWR Method:")
    for idx, row in df_comparison.nsmallest(5, 'rwr_rank').iterrows():
        print(f"  {row['rwr_rank']}. {row['drug']:20s} (score={row['rwr_score']:.6f}, prox_rank={row['proximity_rank']:.0f})")
    
    print("\nProximity Method:")
    for idx, row in df_comparison.nsmallest(5, 'proximity_rank').iterrows():
        print(f"  {row['proximity_rank']}. {row['drug']:20s} (distance={row['distance']:.2f}, rwr_rank={row['rwr_rank']:.0f})")
    
    # Biggest disagreements
    print(f"\n[BIGGEST DISAGREEMENTS]")
    print("Drugs with largest rank differences:")
    for idx, row in df_comparison.nlargest(5, 'rank_diff').iterrows():
        print(f"  • {row['drug']:20s} RWR={row['rwr_rank']:.0f}, Proximity={row['proximity_rank']:.0f}, diff={row['rank_diff']:.0f}")


# ========================================
# Main Execution
# ========================================

def main():
    """Execute bonus RWR analysis."""
    
    print("="*80)
    print(f"BONUS — RANDOM WALK WITH RESTART ({HANDLE})")
    print("="*80)
    
    # Load data
    print("\n[STEP 1] Loading data...")
    G = load_graph(GRAPH_PICKLE)
    print(f"  ✓ Loaded graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    disease_genes = load_disease_genes(DISEASE_GENES_TXT)
    print(f"  ✓ Disease genes: {disease_genes}")
    
    df_proximity = pd.read_csv(DRUG_PRIORITY_CSV)
    print(f"  ✓ Loaded proximity ranking: {len(df_proximity)} drugs")
    
    # Run RWR
    print("\n[STEP 2] Running Random Walk with Restart...")
    rwr_scores = random_walk_with_restart(
        G, disease_genes, 
        restart_prob=RESTART_PROB,
        max_iter=MAX_ITER,
        tol=TOLERANCE
    )
    
    # Rank drugs by RWR
    print("\n[STEP 3] Ranking drugs by RWR scores...")
    df_rwr = rank_drugs_by_rwr(G, rwr_scores)
    df_rwr.to_csv(OUT_RWR_RANKING, index=False)
    print(f"  ✓ Saved RWR ranking to {OUT_RWR_RANKING}")
    
    # Compare rankings
    print("\n[STEP 4] Comparing RWR with Proximity method...")
    df_comparison = compare_rankings(df_rwr, df_proximity)
    df_comparison.to_csv(OUT_COMPARISON, index=False)
    print(f"  ✓ Saved comparison to {OUT_COMPARISON}")
    
    # Analyze comparison
    analyze_comparison(df_comparison)
    
    # Visualize
    print("\n[STEP 5] Creating comparison visualization...")
    visualize_comparison(df_comparison, OUT_COMPARISON_PLOT)
    
    # Summary
    print("\n" + "="*80)
    print("BONUS TASK COMPLETED")
    print("="*80)
    print("\nGenerated files:")
    print(f"  1. {OUT_RWR_RANKING}")
    print(f"  2. {OUT_COMPARISON}")
    print(f"  3. {OUT_COMPARISON_PLOT}")
    print("\n✓ Random Walk with Restart analysis complete!\n")


if __name__ == "__main__":
    main()
