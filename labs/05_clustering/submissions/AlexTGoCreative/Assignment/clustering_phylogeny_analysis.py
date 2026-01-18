#!/usr/bin/env python3
"""
Assignment 5 - Clustering and Phylogenetics Integration
Author: AlexTGoCreative

This script:
1. Loads the multi-FASTA dataset and phylogenetic tree from Lab 4
2. Applies multiple clustering methods (K-means, Hierarchical, DBSCAN)
3. Compares clustering results with phylogenetic tree structure
4. Generates visualizations and analysis report
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import KMeans, AgglomerativeClustering, DBSCAN
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import pandas as pd
from collections import defaultdict

# Set style for better visualizations
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300


def load_sequences(fasta_file):
    """Load sequences from FASTA file"""
    sequences = []
    names = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        names.append(record.id)
    return sequences, names


def sequence_to_kmer_vector(sequence, k=3):
    """Convert sequence to k-mer frequency vector"""
    kmers = {}
    # Generate all possible k-mers
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if 'X' not in kmer and '-' not in kmer:  # Skip ambiguous kmers
            kmers[kmer] = kmers.get(kmer, 0) + 1
    return kmers


def create_kmer_matrix(sequences, k=3):
    """Create k-mer frequency matrix for all sequences"""
    all_kmers = set()
    seq_kmers = []
    
    # Get all k-mers from all sequences
    for seq in sequences:
        kmers = sequence_to_kmer_vector(seq, k)
        seq_kmers.append(kmers)
        all_kmers.update(kmers.keys())
    
    all_kmers = sorted(all_kmers)
    
    # Create matrix
    matrix = np.zeros((len(sequences), len(all_kmers)))
    for i, kmers in enumerate(seq_kmers):
        for j, kmer in enumerate(all_kmers):
            matrix[i, j] = kmers.get(kmer, 0)
    
    # Normalize by sequence length
    for i in range(len(sequences)):
        total = np.sum(matrix[i, :])
        if total > 0:
            matrix[i, :] = matrix[i, :] / total
    
    return matrix, all_kmers


def calculate_pairwise_distances(sequences):
    """Calculate pairwise distances between sequences using Hamming distance"""
    n = len(sequences)
    distances = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            # Align to same length and calculate Hamming distance
            seq1, seq2 = sequences[i], sequences[j]
            max_len = max(len(seq1), len(seq2))
            seq1 = seq1.ljust(max_len, '-')
            seq2 = seq2.ljust(max_len, '-')
            
            dist = sum(c1 != c2 for c1, c2 in zip(seq1, seq2)) / max_len
            distances[i, j] = dist
            distances[j, i] = dist
    
    return distances


def apply_kmeans(data, n_clusters_range=[2, 3, 4]):
    """Apply K-means clustering with different k values"""
    results = {}
    
    for k in n_clusters_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(data)
        
        # Calculate metrics
        if k > 1 and len(set(labels)) > 1:
            silhouette = silhouette_score(data, labels)
            davies_bouldin = davies_bouldin_score(data, labels)
            calinski = calinski_harabasz_score(data, labels)
        else:
            silhouette = davies_bouldin = calinski = 0
        
        results[k] = {
            'labels': labels,
            'silhouette': silhouette,
            'davies_bouldin': davies_bouldin,
            'calinski_harabasz': calinski,
            'centers': kmeans.cluster_centers_
        }
    
    return results


def apply_hierarchical(data, n_clusters_range=[2, 3, 4]):
    """Apply Hierarchical clustering with average linkage"""
    results = {}
    
    for k in n_clusters_range:
        hierarchical = AgglomerativeClustering(n_clusters=k, linkage='average')
        labels = hierarchical.fit_predict(data)
        
        # Calculate metrics
        if k > 1 and len(set(labels)) > 1:
            silhouette = silhouette_score(data, labels)
            davies_bouldin = davies_bouldin_score(data, labels)
            calinski = calinski_harabasz_score(data, labels)
        else:
            silhouette = davies_bouldin = calinski = 0
        
        results[k] = {
            'labels': labels,
            'silhouette': silhouette,
            'davies_bouldin': davies_bouldin,
            'calinski_harabasz': calinski
        }
    
    return results


def apply_dbscan(data, eps_values=[0.3, 0.5, 0.7], min_samples=2):
    """Apply DBSCAN clustering"""
    results = {}
    
    for eps in eps_values:
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        labels = dbscan.fit_predict(data)
        
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        
        # Calculate metrics (only if we have valid clusters)
        if n_clusters > 1 and len(set(labels)) > 1:
            # Filter out noise points for metric calculation
            mask = labels != -1
            if np.sum(mask) > 0:
                silhouette = silhouette_score(data[mask], labels[mask])
            else:
                silhouette = 0
        else:
            silhouette = 0
        
        results[eps] = {
            'labels': labels,
            'n_clusters': n_clusters,
            'n_noise': np.sum(labels == -1),
            'silhouette': silhouette
        }
    
    return results


def extract_phylo_clades(tree):
    """Extract clades from phylogenetic tree"""
    clades = {}
    
    # Get all terminals
    terminals = tree.get_terminals()
    
    # For each pair of terminals, find their common ancestor
    # and group sequences by their closest relatives
    for i, term in enumerate(terminals):
        clades[term.name] = i
    
    return clades


def compare_clustering_phylogeny(labels, names, tree_file):
    """Compare clustering results with phylogenetic tree structure"""
    # Load phylogenetic tree
    tree = Phylo.read(tree_file, "newick")
    
    # Get phylogenetic groupings
    phylo_groups = {}
    
    # Simple approach: for each sequence, find its closest relative
    terminals = tree.get_terminals()
    terminal_dict = {term.name: i for i, term in enumerate(terminals)}
    
    # Create mapping between names and phylo positions
    phylo_labels = np.zeros(len(names))
    for i, name in enumerate(names):
        if name in terminal_dict:
            phylo_labels[i] = terminal_dict[name]
        else:
            # Try to find partial match
            for term_name, idx in terminal_dict.items():
                if name in term_name or term_name in name:
                    phylo_labels[i] = idx
                    break
    
    # Compare with clustering labels
    comparison = pd.DataFrame({
        'Sequence': names,
        'Cluster': labels,
        'Phylo_Position': phylo_labels
    })
    
    return comparison, tree


def plot_dendrogram(data, names, output_file):
    """Plot hierarchical clustering dendrogram"""
    plt.figure(figsize=(12, 6))
    
    # Compute linkage
    Z = linkage(data, method='average')
    
    # Plot dendrogram
    dendrogram(Z, labels=names, leaf_rotation=90)
    plt.title('Hierarchical Clustering Dendrogram (Average Linkage)', fontsize=14, fontweight='bold')
    plt.xlabel('Sequence', fontsize=12)
    plt.ylabel('Distance', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Dendrogram saved to {output_file}")


def plot_pca_clusters(data, labels, names, method_name, output_file):
    """Plot PCA projection with cluster colors"""
    # Apply PCA
    pca = PCA(n_components=2)
    data_pca = pca.fit_transform(data)
    
    plt.figure(figsize=(10, 8))
    
    # Plot each cluster with different color
    unique_labels = set(labels)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_labels)))
    
    for label, color in zip(unique_labels, colors):
        mask = labels == label
        if label == -1:
            # Noise points (for DBSCAN)
            plt.scatter(data_pca[mask, 0], data_pca[mask, 1], 
                       c='gray', marker='x', s=100, alpha=0.5,
                       label='Noise')
        else:
            plt.scatter(data_pca[mask, 0], data_pca[mask, 1], 
                       c=[color], s=150, alpha=0.7,
                       label=f'Cluster {label}', edgecolors='black', linewidth=1)
            
            # Add labels
            for i, name in enumerate(names):
                if mask[i]:
                    plt.annotate(name.split('.')[0].split('|')[-1][:15], 
                               (data_pca[i, 0], data_pca[i, 1]),
                               fontsize=8, alpha=0.8)
    
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)', fontsize=12)
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)', fontsize=12)
    plt.title(f'PCA Projection - {method_name}', fontsize=14, fontweight='bold')
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ PCA plot saved to {output_file}")


def plot_phylogenetic_tree_with_clusters(tree_file, cluster_labels, names, output_file):
    """Plot phylogenetic tree with cluster annotations"""
    tree = Phylo.read(tree_file, "newick")
    
    # Create mapping of names to clusters
    name_to_cluster = dict(zip(names, cluster_labels))
    
    # Color terminals based on clusters
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Draw tree
    Phylo.draw(tree, axes=ax, do_show=False)
    
    # Get terminal positions and color them
    colors = plt.cm.tab10(np.linspace(0, 1, len(set(cluster_labels))))
    
    # Add cluster information to terminal labels
    for terminal in tree.get_terminals():
        for name, cluster in name_to_cluster.items():
            if name in terminal.name or terminal.name in name:
                terminal.name = f"{terminal.name} [C{cluster}]"
                break
    
    # Redraw with updated labels
    plt.close()
    
    fig, ax = plt.subplots(figsize=(14, 10))
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: x.name if x.name else '')
    plt.title('Phylogenetic Tree with Cluster Annotations', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Annotated phylogenetic tree saved to {output_file}")


def generate_report_data(kmeans_results, hierarchical_results, dbscan_results, 
                         comparison_df, names):
    """Generate summary data for the report"""
    report = []
    
    report.append("=" * 70)
    report.append("CLUSTERING AND PHYLOGENETICS ANALYSIS REPORT")
    report.append("=" * 70)
    report.append("")
    
    # Dataset info
    report.append(f"Dataset: {len(names)} p53 protein sequences from multiple species")
    report.append("")
    
    # K-means results
    report.append("--- K-MEANS CLUSTERING ---")
    for k, results in sorted(kmeans_results.items()):
        report.append(f"\nK = {k}:")
        report.append(f"  Silhouette Score: {results['silhouette']:.4f}")
        report.append(f"  Davies-Bouldin Index: {results['davies_bouldin']:.4f}")
        report.append(f"  Calinski-Harabasz Index: {results['calinski_harabasz']:.4f}")
    
    # Hierarchical results
    report.append("\n--- HIERARCHICAL CLUSTERING (Average Linkage) ---")
    for k, results in sorted(hierarchical_results.items()):
        report.append(f"\nK = {k}:")
        report.append(f"  Silhouette Score: {results['silhouette']:.4f}")
        report.append(f"  Davies-Bouldin Index: {results['davies_bouldin']:.4f}")
        report.append(f"  Calinski-Harabasz Index: {results['calinski_harabasz']:.4f}")
    
    # DBSCAN results
    report.append("\n--- DBSCAN CLUSTERING ---")
    for eps, results in sorted(dbscan_results.items()):
        report.append(f"\neps = {eps}:")
        report.append(f"  Number of clusters: {results['n_clusters']}")
        report.append(f"  Noise points: {results['n_noise']}")
        report.append(f"  Silhouette Score: {results['silhouette']:.4f}")
    
    # Best clustering comparison
    report.append("\n" + "=" * 70)
    report.append("COMPARISON WITH PHYLOGENETIC TREE")
    report.append("=" * 70)
    report.append("")
    report.append(comparison_df.to_string())
    
    return "\n".join(report)


def main():
    """Main analysis pipeline"""
    print("\n" + "=" * 70)
    print("Assignment 5: Clustering and Phylogenetics Integration")
    print("=" * 70 + "\n")
    
    # File paths (adjust as needed)
    fasta_file = "/workspaces/bioinf-y4-lab/labs/04_phylogenetics/submissions/AlexTGoCreative/assignment/tp53_multi_species.fasta"
    tree_file = "/workspaces/bioinf-y4-lab/labs/04_phylogenetics/submissions/AlexTGoCreative/assignment/tp53_tree.nwk"
    
    # Load data
    print("📂 Loading sequences from FASTA file...")
    sequences, names = load_sequences(fasta_file)
    print(f"   Loaded {len(sequences)} sequences")
    
    # Create feature matrix using k-mer frequencies
    print("\n🧬 Creating k-mer feature matrix (k=3)...")
    kmer_matrix, kmer_list = create_kmer_matrix(sequences, k=3)
    print(f"   Matrix shape: {kmer_matrix.shape}")
    
    # Apply clustering methods
    print("\n🔍 Applying clustering methods...")
    
    print("   • K-means (k=2, 3, 4)...")
    kmeans_results = apply_kmeans(kmer_matrix, n_clusters_range=[2, 3, 4])
    
    print("   • Hierarchical (k=2, 3, 4)...")
    hierarchical_results = apply_hierarchical(kmer_matrix, n_clusters_range=[2, 3, 4])
    
    print("   • DBSCAN (eps=0.3, 0.5, 0.7)...")
    dbscan_results = apply_dbscan(kmer_matrix, eps_values=[0.3, 0.5, 0.7], min_samples=2)
    
    # Select best clustering (highest silhouette score)
    best_k = max(kmeans_results.keys(), key=lambda k: kmeans_results[k]['silhouette'])
    best_kmeans_labels = kmeans_results[best_k]['labels']
    
    best_h_k = max(hierarchical_results.keys(), key=lambda k: hierarchical_results[k]['silhouette'])
    best_hierarchical_labels = hierarchical_results[best_h_k]['labels']
    
    print(f"\n   Best K-means: k={best_k} (Silhouette: {kmeans_results[best_k]['silhouette']:.4f})")
    print(f"   Best Hierarchical: k={best_h_k} (Silhouette: {hierarchical_results[best_h_k]['silhouette']:.4f})")
    
    # Compare with phylogenetic tree
    print("\n🌳 Comparing clustering with phylogenetic tree...")
    comparison_df, tree = compare_clustering_phylogeny(best_kmeans_labels, names, tree_file)
    
    # Generate visualizations
    print("\n📊 Generating visualizations...")
    
    plot_dendrogram(kmer_matrix, names, "dendrogram.png")
    
    plot_pca_clusters(kmer_matrix, best_kmeans_labels, names, 
                     f"K-means (k={best_k})", "pca_kmeans.png")
    
    plot_pca_clusters(kmer_matrix, best_hierarchical_labels, names, 
                     f"Hierarchical (k={best_h_k})", "pca_hierarchical.png")
    
    # Find best DBSCAN with actual clusters
    valid_dbscan = {eps: res for eps, res in dbscan_results.items() if res['n_clusters'] > 1}
    if valid_dbscan:
        best_eps = max(valid_dbscan.keys(), key=lambda e: valid_dbscan[e]['silhouette'])
        best_dbscan_labels = dbscan_results[best_eps]['labels']
        plot_pca_clusters(kmer_matrix, best_dbscan_labels, names, 
                         f"DBSCAN (eps={best_eps})", "pca_dbscan.png")
    
    plot_phylogenetic_tree_with_clusters(tree_file, best_kmeans_labels, names, 
                                         "phylo_tree_annotated.png")
    
    # Generate report
    print("\n📝 Generating analysis report...")
    report_text = generate_report_data(kmeans_results, hierarchical_results, dbscan_results,
                                      comparison_df, names)
    
    with open("analysis_report.txt", "w") as f:
        f.write(report_text)
    print("   Report saved to analysis_report.txt")
    
    # Save cluster assignments
    comparison_df.to_csv("cluster_assignments.csv", index=False)
    print("   Cluster assignments saved to cluster_assignments.csv")
    
    print("\n" + "=" * 70)
    print("✅ Analysis complete!")
    print("=" * 70 + "\n")
    
    print("Generated files:")
    print("  • dendrogram.png - Hierarchical clustering dendrogram")
    print("  • pca_kmeans.png - PCA plot with K-means clusters")
    print("  • pca_hierarchical.png - PCA plot with hierarchical clusters")
    print("  • pca_dbscan.png - PCA plot with DBSCAN clusters")
    print("  • phylo_tree_annotated.png - Phylogenetic tree with cluster annotations")
    print("  • analysis_report.txt - Detailed analysis results")
    print("  • cluster_assignments.csv - Cluster assignments for each sequence")


if __name__ == "__main__":
    main()
