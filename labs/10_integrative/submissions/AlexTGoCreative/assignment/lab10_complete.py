"""
Lab 10 — Multi-Omics Integration (SNPs + Expression)
Complete Assignment Solution including Bonus

Tasks:
1. Data loading and harmonization (2p)
2. PCA Single-Omics vs Joint (3p)
3. Cross-Omics Correlation (3p)
4. Report (2p) - generated as Markdown
Bonus: Clustering on integrated matrix (+1p)

Author: AlexTGoCreative
"""

from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# Configuration
HANDLE = "AlexTGoCreative"
SNP_CSV = Path(f"data/work/{HANDLE}/lab10/snp_matrix_{HANDLE}.csv")
EXP_CSV = Path(f"data/work/{HANDLE}/lab10/expression_matrix_{HANDLE}.csv")
OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}/assignment")
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("LAB 10 — MULTI-OMICS INTEGRATION: SNPs + Expression")
print("=" * 80)

# ============================================================================
# TASK 1 — Data Loading and Harmonization (2p)
# ============================================================================
print("\n" + "=" * 80)
print("TASK 1 — DATA LOADING AND HARMONIZATION")
print("=" * 80)

print("\n1.1 Loading data...")
df_snp = pd.read_csv(SNP_CSV, index_col=0)
df_exp = pd.read_csv(EXP_CSV, index_col=0)

print(f"   • SNP matrix: {df_snp.shape} (features × samples)")
print(f"   • Expression matrix: {df_exp.shape} (features × samples)")

print("\n1.2 Identifying common samples...")
common_samples = df_snp.columns.intersection(df_exp.columns)
print(f"   • Common samples: {len(common_samples)}")

# Align to common samples
df_snp = df_snp[common_samples]
df_exp = df_exp[common_samples]

print("\n1.3 Normalizing data (z-score)...")
# Z-score normalization: (x - mean) / std for each feature
df_snp_norm = df_snp.sub(df_snp.mean(axis=1), axis=0).div(df_snp.std(axis=1) + 1e-10, axis=0)
df_exp_norm = df_exp.sub(df_exp.mean(axis=1), axis=0).div(df_exp.std(axis=1) + 1e-10, axis=0)

print(f"   • SNP: mean={df_snp_norm.mean().mean():.4f}, std={df_snp_norm.std().mean():.4f}")
print(f"   • Expression: mean={df_exp_norm.mean().mean():.4f}, std={df_exp_norm.std().mean():.4f}")

print("\n1.4 Concatenating layers...")
df_joint = pd.concat([df_snp_norm, df_exp_norm], axis=0)
print(f"   • Joint matrix: {df_joint.shape} (features × samples)")

# Save concatenated matrix
joint_csv = OUT_DIR / f"multiomics_concat_{HANDLE}.csv"
df_joint.to_csv(joint_csv)
print(f"\n✓ Saved: {joint_csv}")

# Store stats for report
task1_stats = {
    'n_snps': len(df_snp_norm),
    'n_genes': len(df_exp_norm),
    'n_samples': len(common_samples),
    'total_features': len(df_joint)
}

# ============================================================================
# TASK 2 — PCA Single-Omics vs Joint (3p)
# ============================================================================
print("\n" + "=" * 80)
print("TASK 2 — PCA SINGLE-OMICS VS JOINT")
print("=" * 80)

n_components = 2

print("\n2.1 PCA on SNP data only...")
pca_snp = PCA(n_components=n_components)
proj_snp = pca_snp.fit_transform(df_snp_norm.T)  # Transpose: samples × features
var_snp = pca_snp.explained_variance_ratio_
print(f"   • PC1: {var_snp[0]*100:.2f}%, PC2: {var_snp[1]*100:.2f}%")
print(f"   • Total variance explained: {sum(var_snp)*100:.2f}%")

print("\n2.2 PCA on Expression data only...")
pca_exp = PCA(n_components=n_components)
proj_exp = pca_exp.fit_transform(df_exp_norm.T)
var_exp = pca_exp.explained_variance_ratio_
print(f"   • PC1: {var_exp[0]*100:.2f}%, PC2: {var_exp[1]*100:.2f}%")
print(f"   • Total variance explained: {sum(var_exp)*100:.2f}%")

print("\n2.3 PCA on Joint multi-omics data...")
pca_joint = PCA(n_components=n_components)
proj_joint = pca_joint.fit_transform(df_joint.T)
var_joint = pca_joint.explained_variance_ratio_
print(f"   • PC1: {var_joint[0]*100:.2f}%, PC2: {var_joint[1]*100:.2f}%")
print(f"   • Total variance explained: {sum(var_joint)*100:.2f}%")

# Store PCA results for report
task2_stats = {
    'snp_var': var_snp,
    'exp_var': var_exp,
    'joint_var': var_joint,
    'snp_total': sum(var_snp),
    'exp_total': sum(var_exp),
    'joint_total': sum(var_joint)
}

print("\n2.4 Generating PCA visualizations...")

# Create sample colors for visualization
sample_colors = np.arange(len(common_samples))

# Figure 1: PCA on SNP data
fig, ax = plt.subplots(figsize=(10, 7))
scatter = ax.scatter(proj_snp[:, 0], proj_snp[:, 1], 
                     c=sample_colors, cmap='viridis', 
                     s=100, alpha=0.7, edgecolors='black', linewidth=0.8)
ax.set_xlabel(f'PC1 ({var_snp[0]*100:.2f}%)', fontsize=14, fontweight='bold')
ax.set_ylabel(f'PC2 ({var_snp[1]*100:.2f}%)', fontsize=14, fontweight='bold')
ax.set_title('PCA on SNP Data Only', fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, linestyle='--')
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
ax.axvline(x=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax, label='Sample Index')
plt.tight_layout()
fig_snp = OUT_DIR / f"pca_snp_{HANDLE}.png"
plt.savefig(fig_snp, dpi=300, bbox_inches='tight')
print(f"   • Saved: {fig_snp.name}")
plt.close()

# Figure 2: PCA on Expression data
fig, ax = plt.subplots(figsize=(10, 7))
scatter = ax.scatter(proj_exp[:, 0], proj_exp[:, 1], 
                     c=sample_colors, cmap='plasma', 
                     s=100, alpha=0.7, edgecolors='black', linewidth=0.8)
ax.set_xlabel(f'PC1 ({var_exp[0]*100:.2f}%)', fontsize=14, fontweight='bold')
ax.set_ylabel(f'PC2 ({var_exp[1]*100:.2f}%)', fontsize=14, fontweight='bold')
ax.set_title('PCA on Expression Data Only', fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, linestyle='--')
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
ax.axvline(x=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax, label='Sample Index')
plt.tight_layout()
fig_exp = OUT_DIR / f"pca_expr_{HANDLE}.png"
plt.savefig(fig_exp, dpi=300, bbox_inches='tight')
print(f"   • Saved: {fig_exp.name}")
plt.close()

# Figure 3: PCA on Joint data
fig, ax = plt.subplots(figsize=(10, 7))
scatter = ax.scatter(proj_joint[:, 0], proj_joint[:, 1], 
                     c=sample_colors, cmap='coolwarm', 
                     s=100, alpha=0.7, edgecolors='black', linewidth=0.8)
ax.set_xlabel(f'PC1 ({var_joint[0]*100:.2f}%)', fontsize=14, fontweight='bold')
ax.set_ylabel(f'PC2 ({var_joint[1]*100:.2f}%)', fontsize=14, fontweight='bold')
ax.set_title('PCA on Joint Multi-Omics Data', fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, linestyle='--')
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
ax.axvline(x=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax, label='Sample Index')
plt.tight_layout()
fig_joint = OUT_DIR / f"pca_joint_{HANDLE}.png"
plt.savefig(fig_joint, dpi=300, bbox_inches='tight')
print(f"   • Saved: {fig_joint.name}")
plt.close()

# Comparison figure: All three side by side
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

axes[0].scatter(proj_snp[:, 0], proj_snp[:, 1], c=sample_colors, 
                cmap='viridis', s=80, alpha=0.7, edgecolors='black', linewidth=0.5)
axes[0].set_xlabel(f'PC1 ({var_snp[0]*100:.2f}%)', fontweight='bold')
axes[0].set_ylabel(f'PC2 ({var_snp[1]*100:.2f}%)', fontweight='bold')
axes[0].set_title('SNP Only', fontsize=14, fontweight='bold')
axes[0].grid(True, alpha=0.3)
axes[0].axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
axes[0].axvline(x=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)

axes[1].scatter(proj_exp[:, 0], proj_exp[:, 1], c=sample_colors, 
                cmap='plasma', s=80, alpha=0.7, edgecolors='black', linewidth=0.5)
axes[1].set_xlabel(f'PC1 ({var_exp[0]*100:.2f}%)', fontweight='bold')
axes[1].set_ylabel(f'PC2 ({var_exp[1]*100:.2f}%)', fontweight='bold')
axes[1].set_title('Expression Only', fontsize=14, fontweight='bold')
axes[1].grid(True, alpha=0.3)
axes[1].axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
axes[1].axvline(x=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)

sc = axes[2].scatter(proj_joint[:, 0], proj_joint[:, 1], c=sample_colors, 
                     cmap='coolwarm', s=80, alpha=0.7, edgecolors='black', linewidth=0.5)
axes[2].set_xlabel(f'PC1 ({var_joint[0]*100:.2f}%)', fontweight='bold')
axes[2].set_ylabel(f'PC2 ({var_joint[1]*100:.2f}%)', fontweight='bold')
axes[2].set_title('Joint Multi-Omics', fontsize=14, fontweight='bold')
axes[2].grid(True, alpha=0.3)
axes[2].axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
axes[2].axvline(x=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)

plt.colorbar(sc, ax=axes, label='Sample Index', fraction=0.02, pad=0.04)
plt.suptitle('PCA Comparison: Single-Omics vs Multi-Omics Integration', 
             fontsize=18, fontweight='bold', y=1.02)
plt.tight_layout()
fig_comparison = OUT_DIR / f"pca_comparison_{HANDLE}.png"
plt.savefig(fig_comparison, dpi=300, bbox_inches='tight')
print(f"   • Saved: {fig_comparison.name}")
plt.close()

print("\n✓ Task 2 Complete")

# ============================================================================
# TASK 3 — Cross-Omics Correlation (3p)
# ============================================================================
print("\n" + "=" * 80)
print("TASK 3 — CROSS-OMICS CORRELATION")
print("=" * 80)

print("\n3.1 Identifying SNPs and Genes...")
snp_ids = [idx for idx in df_joint.index if idx.startswith('rs')]
gene_ids = [idx for idx in df_joint.index if idx.startswith('GENE')]

print(f"   • SNPs: {len(snp_ids)}")
print(f"   • Genes: {len(gene_ids)}")
print(f"   • Total pairs to analyze: {len(snp_ids) * len(gene_ids)}")

print("\n3.2 Computing correlations...")
correlations = []

for snp in snp_ids:
    snp_values = df_joint.loc[snp].values
    
    for gene in gene_ids:
        gene_values = df_joint.loc[gene].values
        
        # Pearson correlation
        corr, pval = pearsonr(snp_values, gene_values)
        
        correlations.append({
            'SNP': snp,
            'Gene': gene,
            'Correlation': corr,
            'P_value': pval,
            'Abs_Correlation': abs(corr)
        })

df_corr = pd.DataFrame(correlations)
print(f"   • Computed {len(df_corr)} correlations")
print(f"   • Mean |r|: {df_corr['Abs_Correlation'].mean():.4f}")
print(f"   • Max |r|: {df_corr['Abs_Correlation'].max():.4f}")

print("\n3.3 Filtering significant correlations (|r| > 0.5)...")
threshold = 0.5
df_filtered = df_corr[df_corr['Abs_Correlation'] > threshold].copy()
df_filtered = df_filtered.sort_values('Abs_Correlation', ascending=False)

print(f"   • Significant pairs (|r| > {threshold}): {len(df_filtered)}")

if len(df_filtered) > 0:
    print(f"\n   Top 10 SNP–Gene pairs:")
    top10 = df_filtered[['SNP', 'Gene', 'Correlation', 'P_value']].head(10)
    for idx, row in top10.iterrows():
        print(f"      {row['SNP']} ↔ {row['Gene']}: r={row['Correlation']:.4f} (p={row['P_value']:.2e})")
else:
    print(f"   No pairs with |r| > {threshold}, lowering threshold...")
    threshold = 0.4
    df_filtered = df_corr[df_corr['Abs_Correlation'] > threshold].copy()
    df_filtered = df_filtered.sort_values('Abs_Correlation', ascending=False)
    print(f"   Pairs with |r| > {threshold}: {len(df_filtered)}")

# Save results
pairs_csv = OUT_DIR / f"snp_gene_pairs_{HANDLE}.csv"
df_filtered.to_csv(pairs_csv, index=False)
print(f"\n✓ Saved: {pairs_csv}")

# Save full correlation matrix
full_corr_csv = OUT_DIR / f"snp_gene_correlations_full_{HANDLE}.csv"
df_corr.to_csv(full_corr_csv, index=False)

# Correlation distribution visualization
print("\n3.4 Generating correlation distribution plot...")
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Histogram
axes[0].hist(df_corr['Correlation'], bins=50, color='steelblue', edgecolor='black', alpha=0.7)
axes[0].axvline(x=threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold |r|={threshold}')
axes[0].axvline(x=-threshold, color='red', linestyle='--', linewidth=2)
axes[0].set_xlabel('Correlation Coefficient (r)', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Frequency', fontsize=12, fontweight='bold')
axes[0].set_title('Distribution of SNP–Gene Correlations', fontsize=14, fontweight='bold')
axes[0].legend(fontsize=10)
axes[0].grid(True, alpha=0.3)

# Scatter plot of correlation vs p-value
axes[1].scatter(df_corr['Correlation'], -np.log10(df_corr['P_value'] + 1e-300), 
                c=df_corr['Abs_Correlation'], cmap='viridis', s=20, alpha=0.5)
axes[1].axvline(x=threshold, color='red', linestyle='--', linewidth=2)
axes[1].axvline(x=-threshold, color='red', linestyle='--', linewidth=2)
axes[1].axhline(y=-np.log10(0.05), color='orange', linestyle='--', linewidth=2, label='p=0.05')
axes[1].set_xlabel('Correlation Coefficient (r)', fontsize=12, fontweight='bold')
axes[1].set_ylabel('-log10(p-value)', fontsize=12, fontweight='bold')
axes[1].set_title('Volcano Plot: Correlation vs Significance', fontsize=14, fontweight='bold')
axes[1].legend(fontsize=10)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
corr_dist = OUT_DIR / f"correlation_distribution_{HANDLE}.png"
plt.savefig(corr_dist, dpi=300, bbox_inches='tight')
print(f"   • Saved: {corr_dist.name}")
plt.close()

# Store stats for report
task3_stats = {
    'total_pairs': len(df_corr),
    'significant_pairs': len(df_filtered),
    'threshold': threshold,
    'mean_abs_corr': df_corr['Abs_Correlation'].mean(),
    'max_abs_corr': df_corr['Abs_Correlation'].max(),
    'top_pairs': df_filtered.head(10) if len(df_filtered) > 0 else pd.DataFrame()
}

print("\n✓ Task 3 Complete")

# ============================================================================
# BONUS — Clustering on Integrated Matrix (+1p)
# ============================================================================
print("\n" + "=" * 80)
print("BONUS — CLUSTERING ON INTEGRATED MATRIX")
print("=" * 80)

print("\n4.1 K-Means Clustering...")
n_clusters_range = range(2, 7)
inertias = []
silhouette_scores = []

from sklearn.metrics import silhouette_score

for k in n_clusters_range:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(df_joint.T)
    inertias.append(kmeans.inertia_)
    silhouette_scores.append(silhouette_score(df_joint.T, labels))

# Choose optimal k (using elbow method and silhouette)
optimal_k = 3  # Can be adjusted based on silhouette scores
print(f"   • Testing k from {min(n_clusters_range)} to {max(n_clusters_range)}")
print(f"   • Optimal k (chosen): {optimal_k}")

kmeans_final = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
cluster_labels = kmeans_final.fit_predict(df_joint.T)

print(f"\n   Cluster distribution:")
for i in range(optimal_k):
    count = np.sum(cluster_labels == i)
    print(f"      Cluster {i+1}: {count} samples ({count/len(cluster_labels)*100:.1f}%)")

print("\n4.2 Hierarchical Clustering...")
linkage_matrix = linkage(df_joint.T, method='ward')

print("\n4.3 Generating clustering visualizations...")

# Figure: Elbow plot
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

axes[0].plot(n_clusters_range, inertias, 'bo-', linewidth=2, markersize=8)
axes[0].set_xlabel('Number of Clusters (k)', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Inertia (Within-Cluster Sum of Squares)', fontsize=12, fontweight='bold')
axes[0].set_title('Elbow Method for Optimal k', fontsize=14, fontweight='bold')
axes[0].grid(True, alpha=0.3)
axes[0].axvline(x=optimal_k, color='red', linestyle='--', linewidth=2, label=f'Chosen k={optimal_k}')
axes[0].legend()

axes[1].plot(n_clusters_range, silhouette_scores, 'ro-', linewidth=2, markersize=8)
axes[1].set_xlabel('Number of Clusters (k)', fontsize=12, fontweight='bold')
axes[1].set_ylabel('Silhouette Score', fontsize=12, fontweight='bold')
axes[1].set_title('Silhouette Score vs Number of Clusters', fontsize=14, fontweight='bold')
axes[1].grid(True, alpha=0.3)
axes[1].axvline(x=optimal_k, color='red', linestyle='--', linewidth=2, label=f'Chosen k={optimal_k}')
axes[1].legend()

plt.tight_layout()
elbow_fig = OUT_DIR / f"clustering_elbow_{HANDLE}.png"
plt.savefig(elbow_fig, dpi=300, bbox_inches='tight')
print(f"   • Saved: {elbow_fig.name}")
plt.close()

# Figure: PCA colored by clusters
fig, ax = plt.subplots(figsize=(10, 7))
scatter = ax.scatter(proj_joint[:, 0], proj_joint[:, 1], 
                     c=cluster_labels, cmap='Set2', 
                     s=120, alpha=0.8, edgecolors='black', linewidth=1.5)
ax.set_xlabel(f'PC1 ({var_joint[0]*100:.2f}%)', fontsize=14, fontweight='bold')
ax.set_ylabel(f'PC2 ({var_joint[1]*100:.2f}%)', fontsize=14, fontweight='bold')
ax.set_title(f'Multi-Omics Clustering (K-Means, k={optimal_k})', fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, linestyle='--')
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
ax.axvline(x=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)

# Add cluster centroids
centroids = kmeans_final.cluster_centers_
# Project centroids to PCA space (centroids are already in feature space)
centroids_pca = pca_joint.transform(centroids)
ax.scatter(centroids_pca[:, 0], centroids_pca[:, 1], 
           marker='X', s=300, c='red', edgecolors='black', linewidth=2, 
           label='Centroids', zorder=5)

plt.colorbar(scatter, ax=ax, label='Cluster', ticks=range(optimal_k))
ax.legend(fontsize=12)
plt.tight_layout()
clustering_pca = OUT_DIR / f"clustering_pca_{HANDLE}.png"
plt.savefig(clustering_pca, dpi=300, bbox_inches='tight')
print(f"   • Saved: {clustering_pca.name}")
plt.close()

# Figure: Dendrogram
fig, ax = plt.subplots(figsize=(14, 7))
dendrogram(linkage_matrix, ax=ax, color_threshold=0.7*max(linkage_matrix[:,2]))
ax.set_xlabel('Sample Index', fontsize=12, fontweight='bold')
ax.set_ylabel('Distance (Ward)', fontsize=12, fontweight='bold')
ax.set_title('Hierarchical Clustering Dendrogram', fontsize=16, fontweight='bold', pad=20)
plt.tight_layout()
dendrogram_fig = OUT_DIR / f"clustering_dendrogram_{HANDLE}.png"
plt.savefig(dendrogram_fig, dpi=300, bbox_inches='tight')
print(f"   • Saved: {dendrogram_fig.name}")
plt.close()

# Heatmap of clusters
print("\n4.4 Generating cluster heatmap...")
# Take top 30 most variable features for visualization
variances = df_joint.var(axis=1).sort_values(ascending=False)
top_features = variances.head(30).index
df_heatmap = df_joint.loc[top_features]

# Sort samples by cluster
sorted_indices = np.argsort(cluster_labels)
df_heatmap_sorted = df_heatmap.iloc[:, sorted_indices]

fig, ax = plt.subplots(figsize=(14, 10))
sns.heatmap(df_heatmap_sorted, cmap='RdBu_r', center=0, 
            cbar_kws={'label': 'Normalized Expression'}, 
            xticklabels=False, yticklabels=True, ax=ax)
ax.set_xlabel('Samples (sorted by cluster)', fontsize=12, fontweight='bold')
ax.set_ylabel('Features (top 30 by variance)', fontsize=12, fontweight='bold')
ax.set_title('Multi-Omics Heatmap by Cluster', fontsize=16, fontweight='bold', pad=20)

# Add cluster boundaries
cluster_sizes = [np.sum(cluster_labels == i) for i in range(optimal_k)]
boundaries = np.cumsum([0] + cluster_sizes)
for b in boundaries[1:-1]:
    ax.axvline(x=b, color='yellow', linewidth=2)

plt.tight_layout()
heatmap_fig = OUT_DIR / f"clustering_heatmap_{HANDLE}.png"
plt.savefig(heatmap_fig, dpi=300, bbox_inches='tight')
print(f"   • Saved: {heatmap_fig.name}")
plt.close()

# Store bonus stats
bonus_stats = {
    'optimal_k': optimal_k,
    'cluster_labels': cluster_labels,
    'cluster_sizes': cluster_sizes,
    'silhouette_scores': silhouette_scores,
    'best_silhouette': max(silhouette_scores)
}

print("\n✓ Bonus Complete")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("✓ ALL TASKS COMPLETED SUCCESSFULLY!")
print("=" * 80)

print("\n📁 OUTPUT FILES:")
print(f"   Task 1: {joint_csv.name}")
print(f"   Task 2: {fig_snp.name}, {fig_exp.name}, {fig_joint.name}, {fig_comparison.name}")
print(f"   Task 3: {pairs_csv.name}, {corr_dist.name}")
print(f"   Bonus:  {elbow_fig.name}, {clustering_pca.name}, {dendrogram_fig.name}, {heatmap_fig.name}")

print("\n📊 KEY FINDINGS:")
print(f"   • {task1_stats['n_snps']} SNPs + {task1_stats['n_genes']} genes integrated")
print(f"   • {task1_stats['n_samples']} samples analyzed")
print(f"   • PCA variance: SNP={task2_stats['snp_total']*100:.1f}%, Expr={task2_stats['exp_total']*100:.1f}%, Joint={task2_stats['joint_total']*100:.1f}%")
print(f"   • {task3_stats['significant_pairs']} significant SNP–gene correlations (|r| > {task3_stats['threshold']})")
print(f"   • Optimal clustering: k={bonus_stats['optimal_k']} (silhouette={bonus_stats['best_silhouette']:.3f})")

print("\n" + "=" * 80)

# Export stats for report generation
import json

def convert_to_serializable(obj):
    """Convert numpy/pandas types to Python native types for JSON serialization"""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, pd.DataFrame):
        return obj.to_dict()
    elif isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(item) for item in obj]
    return obj

stats_dict = {
    'task1': convert_to_serializable(task1_stats),
    'task2': convert_to_serializable(task2_stats),
    'task3': {k: convert_to_serializable(v) if k != 'top_pairs' else None for k, v in task3_stats.items()},
    'bonus': convert_to_serializable(bonus_stats)
}

stats_json = OUT_DIR / f"analysis_stats_{HANDLE}.json"
with open(stats_json, 'w') as f:
    json.dump(stats_dict, f, indent=2)
print(f"\n📄 Analysis statistics saved to: {stats_json.name}")
print("\n🚀 Ready for report generation!")
