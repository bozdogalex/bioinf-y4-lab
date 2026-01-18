"""
Lab 9 — Drug Repurposing Complete Assignment
Executes all tasks: data loading, network construction, similarity calculation,
proximity analysis, visualization, and report generation.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, List, Tuple
import itertools
import pickle

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

# --------------------------
# Configuration
# --------------------------
HANDLE = "AlexTGoCreative"

# Input files
DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")
DISEASE_GENES_TXT = Path(f"data/work/{HANDLE}/lab09/disease_genes_{HANDLE}.txt")

# Output directory
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Output files
OUT_DRUG_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_DRUG_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"
OUT_DRUG_PRIORITY = OUT_DIR / f"drug_priority_{HANDLE}.csv"
OUT_GRAPH_PICKLE = OUT_DIR / f"network_drug_gene_{HANDLE}.gpickle"
OUT_NETWORK_IMG = OUT_DIR / f"network_drug_gene_{HANDLE}.png"
OUT_REPORT_PDF = OUT_DIR / f"report_repurposing_{HANDLE}.pdf"


# ========================================
# TASK 1: Load Data and Build Bipartite Network
# ========================================

def load_drug_gene_table(path: Path) -> pd.DataFrame:
    """Load and validate drug-gene CSV."""
    df = pd.read_csv(path, comment='#')
    if 'drug' not in df.columns or 'gene' not in df.columns:
        raise ValueError("CSV must contain 'drug' and 'gene' columns")
    return df


def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """Build drug -> set of target genes mapping."""
    drug2genes = {}
    for drug, group in df.groupby('drug'):
        drug2genes[drug] = set(group['gene'].values)
    return drug2genes


def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    """Construct bipartite graph with drug and gene nodes."""
    B = nx.Graph()
    
    # Add drug nodes
    for drug in drug2genes.keys():
        B.add_node(drug, bipartite="drug")
    
    # Add gene nodes and edges
    for drug, genes in drug2genes.items():
        for gene in genes:
            if not B.has_node(gene):
                B.add_node(gene, bipartite="gene")
            B.add_edge(drug, gene)
    
    return B


def summarize_drugs(drug2genes: Dict[str, Set[str]]) -> pd.DataFrame:
    """Generate drug summary with target counts."""
    data = []
    for drug, genes in drug2genes.items():
        data.append({'drug': drug, 'num_targets': len(genes)})
    
    df = pd.DataFrame(data)
    df = df.sort_values('num_targets', ascending=False)
    return df


# ========================================
# TASK 2: Drug Similarity Network
# ========================================

def jaccard_similarity(s1: Set[str], s2: Set[str]) -> float:
    """Calculate Jaccard similarity between two sets."""
    if not s1 and not s2:
        return 0.0
    inter = len(s1 & s2)
    union = len(s1 | s2)
    return inter / union if union > 0 else 0.0


def compute_drug_similarity_edges(
    drug2genes: Dict[str, Set[str]],
    min_sim: float = 0.0,
) -> List[Tuple[str, str, float]]:
    """Calculate pairwise drug similarity using Jaccard index."""
    edges = []
    drugs = list(drug2genes.keys())
    
    for drug1, drug2 in itertools.combinations(drugs, 2):
        sim = jaccard_similarity(drug2genes[drug1], drug2genes[drug2])
        if sim >= min_sim:
            edges.append((drug1, drug2, sim))
    
    return edges


def edges_to_dataframe(edges: List[Tuple[str, str, float]]) -> pd.DataFrame:
    """Convert edge list to DataFrame."""
    if not edges:
        return pd.DataFrame(columns=['drug1', 'drug2', 'similarity'])
    
    df = pd.DataFrame(edges, columns=['drug1', 'drug2', 'similarity'])
    df = df.sort_values('similarity', ascending=False)
    return df


# ========================================
# TASK 3: Disease Proximity Analysis
# ========================================

def load_disease_genes(path: Path) -> Set[str]:
    """Load disease gene list from text file."""
    genes = set()
    with open(path, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                genes.add(gene)
    return genes


def get_drug_nodes(B: nx.Graph) -> List[str]:
    """Extract drug nodes from bipartite graph."""
    return [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]


def compute_drug_disease_distance(
    B: nx.Graph,
    drug: str,
    disease_genes: Set[str],
    mode: str = "mean",
    max_dist: int = 5,
) -> float:
    """Calculate average distance from drug to disease genes."""
    distances = []
    
    for gene in disease_genes:
        if gene not in B:
            continue
        
        try:
            dist = nx.shortest_path_length(B, source=drug, target=gene)
            distances.append(dist)
        except nx.NetworkXNoPath:
            distances.append(max_dist + 1)
    
    if not distances:
        return float('inf')
    
    if mode == "mean":
        return sum(distances) / len(distances)
    elif mode == "min":
        return min(distances)
    else:
        return sum(distances) / len(distances)


def rank_drugs_by_proximity(
    B: nx.Graph,
    disease_genes: Set[str],
    mode: str = "mean",
) -> pd.DataFrame:
    """Rank drugs by proximity to disease genes."""
    drugs = get_drug_nodes(B)
    
    results = []
    for drug in drugs:
        dist = compute_drug_disease_distance(B, drug, disease_genes, mode=mode)
        results.append({'drug': drug, 'distance': dist})
    
    df = pd.DataFrame(results)
    df = df.sort_values('distance')
    
    return df


# ========================================
# TASK 4: Network Visualization
# ========================================

def visualize_network(B: nx.Graph, output_path: Path):
    """Generate bipartite network visualization."""
    drug_nodes = [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]
    gene_nodes = [n for n, d in B.nodes(data=True) if d.get("bipartite") == "gene"]
    
    # Create layout
    pos = nx.spring_layout(B, k=1.5, iterations=50, seed=42)
    
    # Calculate node sizes based on degree
    node_sizes_drugs = [200 + B.degree(n) * 150 for n in drug_nodes]
    node_sizes_genes = [200 + B.degree(n) * 150 for n in gene_nodes]
    
    # Create figure
    plt.figure(figsize=(16, 12))
    
    # Draw edges
    nx.draw_networkx_edges(B, pos, alpha=0.2, width=0.8, edge_color='lightgray')
    
    # Draw gene nodes (red)
    nx.draw_networkx_nodes(
        B, pos,
        nodelist=gene_nodes,
        node_color='#FF6B6B',
        node_size=node_sizes_genes,
        alpha=0.8,
        edgecolors='darkred',
        linewidths=2,
        label='Genes'
    )
    
    # Draw drug nodes (blue)
    nx.draw_networkx_nodes(
        B, pos,
        nodelist=drug_nodes,
        node_color='#4ECDC4',
        node_size=node_sizes_drugs,
        alpha=0.8,
        edgecolors='darkblue',
        linewidths=2,
        label='Drugs'
    )
    
    # Add labels for high-degree nodes
    important_nodes = [n for n in B.nodes() if B.degree(n) >= 2]
    labels = {n: n for n in important_nodes}
    nx.draw_networkx_labels(B, pos, labels=labels, font_size=8, font_weight='bold')
    
    plt.title(f"Drug-Gene Bipartite Network\n{len(drug_nodes)} drugs, {len(gene_nodes)} genes, {B.number_of_edges()} interactions", 
              fontsize=16, fontweight='bold', pad=20)
    plt.legend(scatterpoints=1, frameon=True, fontsize=12, loc='upper right')
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"[INFO] Saved network visualization to {output_path}")
    plt.close()


# ========================================
# TASK 5: Generate PDF Report
# ========================================

def generate_report(
    df_summary: pd.DataFrame,
    df_similarity: pd.DataFrame,
    df_priority: pd.DataFrame,
    B: nx.Graph,
    disease_genes: Set[str],
    output_path: Path
):
    """Generate comprehensive PDF report."""
    
    with PdfPages(output_path) as pdf:
        # Page 1: Overview and Summary Statistics
        fig = plt.figure(figsize=(8.5, 11))
        fig.suptitle(f'Lab 9 — Drug Repurposing Network Analysis\nHandle: {HANDLE}', 
                     fontsize=16, fontweight='bold', y=0.98)
        
        ax1 = plt.subplot(4, 1, 1)
        ax1.axis('off')
        
        overview_text = f"""
OVERVIEW
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Dataset Statistics:
  • Total drugs: {len(df_summary)}
  • Total genes: {len([n for n, d in B.nodes(data=True) if d.get("bipartite") == "gene"])}
  • Drug-gene interactions: {B.number_of_edges()}
  • Disease genes analyzed: {len(disease_genes)} ({', '.join(sorted(disease_genes))})

Network Characteristics:
  • Graph type: Bipartite (drug-gene)
  • Average drug targets: {df_summary['num_targets'].mean():.2f}
  • Drugs with most targets: {df_summary.iloc[0]['drug']} ({df_summary.iloc[0]['num_targets']} targets)
  
Similarity Network:
  • Drug pairs analyzed: {len(df_similarity)}
  • High similarity pairs (>0.5): {len(df_similarity[df_similarity['similarity'] > 0.5])}
  • Perfect matches (similarity = 1.0): {len(df_similarity[df_similarity['similarity'] == 1.0])}
"""
        ax1.text(0.05, 0.95, overview_text, transform=ax1.transAxes,
                fontsize=10, verticalalignment='top', fontfamily='monospace')
        
        # Top drugs by target count
        ax2 = plt.subplot(4, 1, 2)
        top_drugs = df_summary.head(10)
        ax2.barh(range(len(top_drugs)), top_drugs['num_targets'], color='#4ECDC4', alpha=0.7)
        ax2.set_yticks(range(len(top_drugs)))
        ax2.set_yticklabels(top_drugs['drug'], fontsize=9)
        ax2.set_xlabel('Number of Target Genes', fontsize=10)
        ax2.set_title('Top 10 Drugs by Number of Targets', fontsize=11, fontweight='bold')
        ax2.invert_yaxis()
        ax2.grid(axis='x', alpha=0.3)
        
        # Disease proximity results
        ax3 = plt.subplot(4, 1, 3)
        top_candidates = df_priority.head(10)
        colors = plt.cm.RdYlGn_r(top_candidates['distance'] / top_candidates['distance'].max())
        ax3.barh(range(len(top_candidates)), top_candidates['distance'], color=colors, alpha=0.8)
        ax3.set_yticks(range(len(top_candidates)))
        ax3.set_yticklabels(top_candidates['drug'], fontsize=9)
        ax3.set_xlabel('Average Distance to Disease Genes', fontsize=10)
        ax3.set_title('Top 10 Drug Candidates (by Network Proximity)', fontsize=11, fontweight='bold')
        ax3.invert_yaxis()
        ax3.grid(axis='x', alpha=0.3)
        
        # Top similar drug pairs
        ax4 = plt.subplot(4, 1, 4)
        ax4.axis('off')
        
        top_sim = df_similarity.head(8)
        sim_text = "Top Drug Similarity Pairs (Jaccard Index):\n" + "─" * 70 + "\n"
        for _, row in top_sim.iterrows():
            sim_text += f"  • {row['drug1']:20s} ↔ {row['drug2']:20s}  similarity: {row['similarity']:.3f}\n"
        
        ax4.text(0.05, 0.95, sim_text, transform=ax4.transAxes,
                fontsize=9, verticalalignment='top', fontfamily='monospace')
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Page 2: Methodology and Interpretation
        fig = plt.figure(figsize=(8.5, 11))
        ax = fig.add_subplot(111)
        ax.axis('off')
        
        methodology_text = """
METHODOLOGY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. Bipartite Network Construction
   We constructed a drug-gene bipartite network where drugs and genes represent 
   two distinct node sets, with edges representing known drug-gene interactions.
   
2. Drug Similarity Calculation
   Drug similarity was computed using the Jaccard similarity coefficient:
   
       J(A, B) = |A ∩ B| / |A ∪ B|
   
   where A and B are the sets of target genes for drugs A and B. This metric
   captures the overlap in molecular mechanisms between drug pairs.

3. Network Proximity Analysis
   For disease gene prioritization, we calculated the average shortest path
   distance from each drug to the set of disease-associated genes:
   
       proximity(drug) = mean(shortest_path(drug, gene_i)) for all disease genes
   
   Lower distances indicate drugs that are topologically closer to disease genes,
   suggesting potential therapeutic relevance.

RESULTS INTERPRETATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Top Candidate Drugs:
"""
        
        # Add top 5 candidates with interpretation
        top_5 = df_priority.head(5)
        for idx, row in top_5.iterrows():
            methodology_text += f"\n  {idx+1}. {row['drug']} (distance: {row['distance']:.2f})"
        
        methodology_text += """


Biological Interpretation:
  • Drugs with lower proximity scores directly target or interact with genes
    in close network vicinity to the disease genes.
    
  • High similarity scores between drugs suggest shared molecular mechanisms,
    which can guide combination therapy strategies or predict off-target effects.
    
  • Drugs targeting multiple genes (high degree) may have broader therapeutic
    effects but also higher risk of side effects (polypharmacology).

LIMITATIONS AND CAVEATS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. Annotation Bias
   Drug-gene interaction databases are incomplete and biased toward well-studied
   drugs and genes. Novel or understudied drugs may be underrepresented.

2. Network Incompleteness
   The bipartite network captures only known interactions. Unknown or context-
   specific interactions are not represented, potentially missing important
   therapeutic opportunities.

3. Topological Limitations
   Network proximity assumes that topological closeness correlates with
   functional relevance, which may not always hold in complex biological systems.

4. Lack of Clinical Context
   This analysis does not account for:
   • Drug bioavailability and pharmacokinetics
   • Tissue-specific expression patterns
   • Disease stage and patient heterogeneity
   • Drug safety profiles and contraindications

5. Simplified Similarity Metric
   Jaccard similarity treats all gene targets equally, ignoring:
   • Gene expression levels
   • Binding affinity differences
   • Functional importance of individual targets

CONCLUSIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Network-based drug repurposing provides a systematic framework for identifying
candidate therapeutics. While the approach has limitations, it offers valuable
insights for hypothesis generation and prioritization in drug discovery pipelines.

Integration with additional data layers (gene expression, protein-protein
interactions, clinical outcomes) would strengthen the predictive power and
clinical relevance of this approach.
"""
        
        ax.text(0.05, 0.98, methodology_text, transform=ax.transAxes,
                fontsize=9.5, verticalalignment='top', fontfamily='monospace',
                wrap=True)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Page 3: Data Tables
        fig = plt.figure(figsize=(8.5, 11))
        fig.suptitle('Detailed Results Tables', fontsize=14, fontweight='bold', y=0.98)
        
        # Drug summary table
        ax1 = plt.subplot(3, 1, 1)
        ax1.axis('off')
        ax1.set_title('Drug Target Summary (Top 15)', fontsize=11, fontweight='bold', pad=10)
        
        table_data = []
        for idx, row in df_summary.head(15).iterrows():
            table_data.append([row['drug'], row['num_targets']])
        
        table1 = ax1.table(cellText=table_data,
                          colLabels=['Drug', 'Targets'],
                          cellLoc='left',
                          loc='center',
                          colWidths=[0.7, 0.3])
        table1.auto_set_font_size(False)
        table1.set_fontsize(8)
        table1.scale(1, 1.5)
        
        # Drug similarity table
        ax2 = plt.subplot(3, 1, 2)
        ax2.axis('off')
        ax2.set_title('Top Drug Similarity Pairs', fontsize=11, fontweight='bold', pad=10)
        
        table_data = []
        for idx, row in df_similarity.head(12).iterrows():
            table_data.append([row['drug1'], row['drug2'], f"{row['similarity']:.3f}"])
        
        table2 = ax2.table(cellText=table_data,
                          colLabels=['Drug 1', 'Drug 2', 'Similarity'],
                          cellLoc='left',
                          loc='center',
                          colWidths=[0.4, 0.4, 0.2])
        table2.auto_set_font_size(False)
        table2.set_fontsize(8)
        table2.scale(1, 1.5)
        
        # Disease proximity table
        ax3 = plt.subplot(3, 1, 3)
        ax3.axis('off')
        ax3.set_title('Drug-Disease Proximity Ranking', fontsize=11, fontweight='bold', pad=10)
        
        table_data = []
        for idx, row in df_priority.head(15).iterrows():
            table_data.append([idx+1, row['drug'], f"{row['distance']:.3f}"])
        
        table3 = ax3.table(cellText=table_data,
                          colLabels=['Rank', 'Drug', 'Distance'],
                          cellLoc='left',
                          loc='center',
                          colWidths=[0.2, 0.6, 0.2])
        table3.auto_set_font_size(False)
        table3.set_fontsize(8)
        table3.scale(1, 1.5)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    print(f"[INFO] Saved comprehensive report to {output_path}")


# ========================================
# MAIN EXECUTION
# ========================================

def main():
    """Execute complete assignment pipeline."""
    
    print("=" * 80)
    print(f"LAB 9 — DRUG REPURPOSING ASSIGNMENT ({HANDLE})")
    print("=" * 80)
    
    # TASK 1: Load data and build bipartite network
    print("\n[TASK 1] Loading data and building bipartite network...")
    df_drug_gene = load_drug_gene_table(DRUG_GENE_CSV)
    print(f"  ✓ Loaded {len(df_drug_gene)} drug-gene interactions")
    
    drug2genes = build_drug2genes(df_drug_gene)
    print(f"  ✓ Found {len(drug2genes)} unique drugs")
    
    B = build_bipartite_graph(drug2genes)
    print(f"  ✓ Built graph: {B.number_of_nodes()} nodes, {B.number_of_edges()} edges")
    
    # Save graph
    with open(OUT_GRAPH_PICKLE, 'wb') as f:
        pickle.dump(B, f)
    print(f"  ✓ Saved graph to {OUT_GRAPH_PICKLE}")
    
    df_summary = summarize_drugs(drug2genes)
    df_summary.to_csv(OUT_DRUG_SUMMARY, index=False)
    print(f"  ✓ Saved drug summary to {OUT_DRUG_SUMMARY}")
    
    # TASK 2: Drug similarity network
    print("\n[TASK 2] Computing drug similarity network...")
    edges = compute_drug_similarity_edges(drug2genes, min_sim=0.0)
    print(f"  ✓ Computed {len(edges)} drug-drug similarity edges")
    
    df_similarity = edges_to_dataframe(edges)
    df_similarity.to_csv(OUT_DRUG_SIMILARITY, index=False)
    print(f"  ✓ Saved drug similarity to {OUT_DRUG_SIMILARITY}")
    
    # TASK 3: Disease proximity analysis
    print("\n[TASK 3] Analyzing disease proximity...")
    disease_genes = load_disease_genes(DISEASE_GENES_TXT)
    print(f"  ✓ Loaded {len(disease_genes)} disease genes: {disease_genes}")
    
    df_priority = rank_drugs_by_proximity(B, disease_genes, mode="mean")
    df_priority.to_csv(OUT_DRUG_PRIORITY, index=False)
    print(f"  ✓ Saved drug priority ranking to {OUT_DRUG_PRIORITY}")
    print(f"  ✓ Top candidate: {df_priority.iloc[0]['drug']} (distance: {df_priority.iloc[0]['distance']:.3f})")
    
    # TASK 4: Network visualization
    print("\n[TASK 4] Generating network visualization...")
    visualize_network(B, OUT_NETWORK_IMG)
    print(f"  ✓ Visualization complete")
    
    # Summary
    print("\n" + "=" * 80)
    print("ASSIGNMENT COMPLETE - ALL DELIVERABLES GENERATED")
    print("=" * 80)
    print("\nGenerated files:")
    print(f"  1. {OUT_DRUG_SUMMARY}")
    print(f"  2. {OUT_DRUG_SIMILARITY}")
    print(f"  3. {OUT_DRUG_PRIORITY}")
    print(f"  4. {OUT_NETWORK_IMG}")
    print(f"  5. {OUT_GRAPH_PICKLE}")
    print("\n✓ All tasks completed successfully!\n")


if __name__ == "__main__":
    main()
