"""
Exercise 10 — PCA Single-Omics vs Joint

TODO:
- încărcați SNP și Expression
- normalizați fiecare strat (z-score)
- rulați PCA pe:
    1) strat SNP
    2) strat Expression
    3) strat Joint (concat)
- generați 3 figuri PNG
- comparați vizual distribuția probelor
"""

from pathlib import Path
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np

HANDLE = "AlexTGoCreative"

SNP_CSV = Path(f"data/work/{HANDLE}/lab10/snp_matrix_{HANDLE}.csv")
EXP_CSV = Path(f"data/work/{HANDLE}/lab10/expression_matrix_{HANDLE}.csv")

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("Exercise 10.1 — PCA Single-Omics vs Joint")
print("=" * 60)

# Load data
print("\n1. Loading data...")
df_snp = pd.read_csv(SNP_CSV, index_col=0)
df_exp = pd.read_csv(EXP_CSV, index_col=0)

print(f"   SNP matrix: {df_snp.shape} (features × samples)")
print(f"   Expression matrix: {df_exp.shape} (features × samples)")

# Find common samples
common_samples = df_snp.columns.intersection(df_exp.columns)
print(f"\n2. Aligning samples...")
print(f"   Common samples: {len(common_samples)}")

df_snp = df_snp[common_samples]
df_exp = df_exp[common_samples]

# Normalize each layer using z-score normalization
print("\n3. Normalizing data (z-score)...")
# Normalize across samples (row-wise) for each feature
df_snp_norm = df_snp.sub(df_snp.mean(axis=1), axis=0).div(df_snp.std(axis=1) + 1e-10, axis=0)
df_exp_norm = df_exp.sub(df_exp.mean(axis=1), axis=0).div(df_exp.std(axis=1) + 1e-10, axis=0)

print(f"   SNP normalized: mean={df_snp_norm.mean().mean():.4f}, std={df_snp_norm.std().mean():.4f}")
print(f"   Expression normalized: mean={df_exp_norm.mean().mean():.4f}, std={df_exp_norm.std().mean():.4f}")

# Concatenate for joint analysis
df_joint = pd.concat([df_snp_norm, df_exp_norm], axis=0)
print(f"\n4. Creating joint matrix...")
print(f"   Joint matrix: {df_joint.shape} (features × samples)")

# Save joint matrix
joint_csv = OUT_DIR / f"multiomics_concat_{HANDLE}.csv"
df_joint.to_csv(joint_csv)
print(f"   Saved: {joint_csv}")

# Run PCA on three scenarios
print("\n5. Running PCA analysis...")
n_components = 2

# PCA 1: SNP only
pca_snp = PCA(n_components=n_components)
proj_snp = pca_snp.fit_transform(df_snp_norm.T)  # Transpose: samples × features
var_snp = pca_snp.explained_variance_ratio_

print(f"   PCA on SNP: PC1={var_snp[0]*100:.2f}%, PC2={var_snp[1]*100:.2f}%")

# PCA 2: Expression only
pca_exp = PCA(n_components=n_components)
proj_exp = pca_exp.fit_transform(df_exp_norm.T)
var_exp = pca_exp.explained_variance_ratio_

print(f"   PCA on Expression: PC1={var_exp[0]*100:.2f}%, PC2={var_exp[1]*100:.2f}%")

# PCA 3: Joint (concatenated)
pca_joint = PCA(n_components=n_components)
proj_joint = pca_joint.fit_transform(df_joint.T)
var_joint = pca_joint.explained_variance_ratio_

print(f"   PCA on Joint: PC1={var_joint[0]*100:.2f}%, PC2={var_joint[1]*100:.2f}%")

# Generate visualizations
print("\n6. Generating visualizations...")

# Figure 1: PCA on SNP data
fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(proj_snp[:, 0], proj_snp[:, 1], 
                     c=range(len(proj_snp)), cmap='viridis', 
                     s=80, alpha=0.7, edgecolors='black', linewidth=0.5)
ax.set_xlabel(f'PC1 ({var_snp[0]*100:.2f}%)', fontsize=12)
ax.set_ylabel(f'PC2 ({var_snp[1]*100:.2f}%)', fontsize=12)
ax.set_title('PCA on SNP Data Only', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
plt.colorbar(scatter, ax=ax, label='Sample Index')
plt.tight_layout()
fig_snp = OUT_DIR / f"pca_snp_{HANDLE}.png"
plt.savefig(fig_snp, dpi=300)
print(f"   Saved: {fig_snp}")
plt.close()

# Figure 2: PCA on Expression data
fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(proj_exp[:, 0], proj_exp[:, 1], 
                     c=range(len(proj_exp)), cmap='plasma', 
                     s=80, alpha=0.7, edgecolors='black', linewidth=0.5)
ax.set_xlabel(f'PC1 ({var_exp[0]*100:.2f}%)', fontsize=12)
ax.set_ylabel(f'PC2 ({var_exp[1]*100:.2f}%)', fontsize=12)
ax.set_title('PCA on Expression Data Only', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
plt.colorbar(scatter, ax=ax, label='Sample Index')
plt.tight_layout()
fig_exp = OUT_DIR / f"pca_expression_{HANDLE}.png"
plt.savefig(fig_exp, dpi=300)
print(f"   Saved: {fig_exp}")
plt.close()

# Figure 3: PCA on Joint data
fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(proj_joint[:, 0], proj_joint[:, 1], 
                     c=range(len(proj_joint)), cmap='coolwarm', 
                     s=80, alpha=0.7, edgecolors='black', linewidth=0.5)
ax.set_xlabel(f'PC1 ({var_joint[0]*100:.2f}%)', fontsize=12)
ax.set_ylabel(f'PC2 ({var_joint[1]*100:.2f}%)', fontsize=12)
ax.set_title('PCA on Joint Multi-Omics Data', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
plt.colorbar(scatter, ax=ax, label='Sample Index')
plt.tight_layout()
fig_joint = OUT_DIR / f"pca_joint_{HANDLE}.png"
plt.savefig(fig_joint, dpi=300)
print(f"   Saved: {fig_joint}")
plt.close()

# Comparison figure: All three in one
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

axes[0].scatter(proj_snp[:, 0], proj_snp[:, 1], c=range(len(proj_snp)), 
                cmap='viridis', s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
axes[0].set_xlabel(f'PC1 ({var_snp[0]*100:.2f}%)')
axes[0].set_ylabel(f'PC2 ({var_snp[1]*100:.2f}%)')
axes[0].set_title('SNP Only')
axes[0].grid(True, alpha=0.3)

axes[1].scatter(proj_exp[:, 0], proj_exp[:, 1], c=range(len(proj_exp)), 
                cmap='plasma', s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
axes[1].set_xlabel(f'PC1 ({var_exp[0]*100:.2f}%)')
axes[1].set_ylabel(f'PC2 ({var_exp[1]*100:.2f}%)')
axes[1].set_title('Expression Only')
axes[1].grid(True, alpha=0.3)

sc = axes[2].scatter(proj_joint[:, 0], proj_joint[:, 1], c=range(len(proj_joint)), 
                     cmap='coolwarm', s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
axes[2].set_xlabel(f'PC1 ({var_joint[0]*100:.2f}%)')
axes[2].set_ylabel(f'PC2 ({var_joint[1]*100:.2f}%)')
axes[2].set_title('Joint Multi-Omics')
axes[2].grid(True, alpha=0.3)

plt.colorbar(sc, ax=axes, label='Sample Index', fraction=0.02)
plt.suptitle('PCA Comparison: Single-Omics vs Joint', fontsize=16, fontweight='bold')
plt.tight_layout()
fig_comparison = OUT_DIR / f"pca_comparison_{HANDLE}.png"
plt.savefig(fig_comparison, dpi=300)
print(f"   Saved: {fig_comparison}")
plt.close()

print("\n" + "=" * 60)
print("✓ Exercise 10.1 Complete!")
print("=" * 60)
print(f"\nOutputs:")
print(f"  - {joint_csv}")
print(f"  - {fig_snp}")
print(f"  - {fig_exp}")
print(f"  - {fig_joint}")
print(f"  - {fig_comparison}")
print("\nInterpretation:")
print(f"  • SNP PCA explains {var_snp[0]*100:.1f}% + {var_snp[1]*100:.1f}% = {sum(var_snp)*100:.1f}% variance")
print(f"  • Expression PCA explains {var_exp[0]*100:.1f}% + {var_exp[1]*100:.1f}% = {sum(var_exp)*100:.1f}% variance")
print(f"  • Joint PCA explains {var_joint[0]*100:.1f}% + {var_joint[1]*100:.1f}% = {sum(var_joint)*100:.1f}% variance")
print(f"  • Joint analysis captures complementary information from both omics layers")
