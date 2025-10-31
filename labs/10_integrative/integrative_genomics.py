import pandas as pd
import numpy as np
import requests
import io
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
BASE_URL = "https://raw.githubusercontent.com/<your_username>/lab9-datasets/main/"
SNP_FILE = "../../data/sample/snp_data.csv"
EXPR_FILE = "../../data/sample/expression_data.csv"
PROT_FILE = "../../data/sample/proteomics_data.csv"
PHENO_FILE = "../../data/sample/phenotypes.csv"

# ------------------------------------------------------------
# Helper function to fetch CSV from GitHub
# ------------------------------------------------------------
def fetch_csv(url):
    r = requests.get(url)
    r.raise_for_status()
    data = r.text
    df = pd.read_csv(io.StringIO(data), index_col=0)
    return df

# ------------------------------------------------------------
# Step 1: Load Data
# ------------------------------------------------------------
snp_url = BASE_URL + SNP_FILE
#expand from 1000 genomics project (http://www.internationalgenome.org/data) or TCGA (https://portal.gdc.cancer.gov/)
expr_url = BASE_URL + EXPR_FILE
#expand from GEO (Gene Expression Omnibus) (https://www.ncbi.nlm.nih.gov/geo/) or GTEx ( https://gtexportal.org/home/datasets)
prot_url = BASE_URL + PROT_FILE
#expand from CPTAC (Clinical Proteomic Tumor Analysis Consortium) Website: https://proteomics.cancer.gov/programs/cptac)
pheno_url = BASE_URL + PHENO_FILE
#expand from TCGA Clinical Data: (https://portal.gdc.cancer.gov/);  For instance, choose a cohort of breast cancer patients and extract their clinical attributes such as survival, drug response, or tumor stage.


snp_df = fetch_csv(snp_url)      # Samples x SNPs
expr_df = fetch_csv(expr_url)    # Genes x Samples
prot_df = fetch_csv(prot_url)    # Proteins x Samples
pheno_df = fetch_csv(pheno_url)  # Samples x Phenotype

# ------------------------------------------------------------
# Step 2: Data Preprocessing
# ------------------------------------------------------------
# Log2 transform for RNA-Seq expression (assuming raw counts)
expr_df = expr_df.apply(lambda x: np.log2(x + 1))

# Z-score normalization for proteomics data
prot_df = pd.DataFrame(StandardScaler().fit_transform(prot_df.T).T,
                       index=prot_df.index, columns=prot_df.columns)

# Align all datasets by common samples
common_samples = snp_df.index.intersection(expr_df.columns).intersection(prot_df.columns).intersection(pheno_df.index)
snp_df = snp_df.loc[common_samples]
expr_df = expr_df[common_samples]
prot_df = prot_df[common_samples]
y = pheno_df.loc[common_samples, 'phenotype']

# ------------------------------------------------------------
# Step 3: Feature Selection
# ------------------------------------------------------------
# Select top variable genes
gene_variances = expr_df.var(axis=1)
top_genes = gene_variances.sort_values(ascending=False).head(200).index
expr_top = expr_df.loc[top_genes]

# Select top variable proteins
protein_variances = prot_df.var(axis=1)
top_proteins = protein_variances.sort_values(ascending=False).head(100).index
prot_top = prot_df.loc[top_proteins]

# ------------------------------------------------------------
# Step 4: Integration
# ------------------------------------------------------------
# Integration strategy: Combine SNP, expression, and proteomics data horizontally
# After feature selection, we have smaller sets of genes and proteins
X = pd.concat([expr_top.T, prot_top.T, snp_df], axis=1, join='inner')
# Now X has samples as rows, and features from expression, proteomics, and SNPs as columns

# Convert phenotype to binary or numeric if needed
# For example, if phenotype is "responder"/"non_responder", we can map them:
y_mapped = y.map({'responder': 1, 'non_responder': 0})

# ------------------------------------------------------------
# Step 5: Dimensionality Reduction & Clustering
# ------------------------------------------------------------
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

pca = PCA(n_components=5)
X_pca = pca.fit_transform(X_scaled)

kmeans = KMeans(n_clusters=2, random_state=42)
clusters = kmeans.fit_predict(X_pca)

# ------------------------------------------------------------
# Step 6: Machine Learning Prediction
# ------------------------------------------------------------
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_pca, y_mapped)
y_pred = rf.predict(X_pca)
accuracy = (y_pred == y_mapped).mean()
print(f"Training accuracy: {accuracy:.2f}")

# ------------------------------------------------------------
# Step 7: Visualization
# ------------------------------------------------------------
plt.figure(figsize=(8,6))
plt.scatter(X_pca[:,0], X_pca[:,1], c=clusters, cmap='viridis', edgecolor='k')
plt.title("PCA of Integrated Data (SNP + Expression + Proteomics)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.colorbar(label="Cluster")
plt.tight_layout()
plt.savefig("pca_clusters_realistic_data.png")
plt.show()

# ------------------------------------------------------------

