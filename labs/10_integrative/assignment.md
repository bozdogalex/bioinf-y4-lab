# Lab 10 — Multi-Omics Integration (SNPs + Expression)

## Obiective
- Construirea unui pipeline integrativ pe două straturi omice.  
- Compararea analizei single-omics versus integrative.  
- Identificarea candidaților SNP–gene relevanți biologic.  
- Vizualizare și interpretare.

---

## Task 1 — Încărcarea și armonizarea datelor (2p)
- încărcați `snp_matrix_<handle>.csv` și `expression_matrix_<handle>.csv`,  
- identificați probele comune,  
- normalizați fiecare strat (z-score),  
- concatenați straturile într-o singură matrice,  
- salvați: `multiomics_concat_<handle>.csv`.

---

## Task 2 — PCA Single-Omics vs Joint (3p)
- rulați PCA doar pe SNPs,  
- rulați PCA doar pe expresie,  
- rulați PCA pe matricea integrată,  
- vizualizați probele colorate după un factor (ex. subtype).  
Fișiere:  
- `pca_snp_<handle>.png`  
- `pca_expr_<handle>.png`  
- `pca_joint_<handle>.png`

---

## Task 3 — Cross-Omics Correlation (3p)
- calculați corelații între fiecare SNP și fiecare genă,  
- păstrați doar |r| > 0.5,  
- exportați: `snp_gene_pairs_<handle>.csv`.

---

## Task 4 — Raport (2p)
Raport PDF (max 3 pagini):
- cum diferă PCA single-omics vs multi-omics,  
- care perechi SNP–gene ies relevante și de ce,  
- potențiale implicații în oncologie sau farmacogenomică,  
- limitări.

---

## Bonus (+1p)
Aplicați clustering (KMeans sau Hierarchical) pe matricea integrată și comparați cu subtipurile clinice.

