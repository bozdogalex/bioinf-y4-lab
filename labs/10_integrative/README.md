# Săptămâna 10 — Integrative Genomics Approaches  

> **Predare 2026–2027:** toate soluțiile și PR-urile de evaluare rămân în repository-ul vostru **privat**, partajat cu `bozdogalex`. Trimiteți linkul PR și commit SHA în LMS. [Ghid de predare](../../docs/git-workflow.md).
## Multi-Omics Fusion for Precision Medicine

## Scopuri
- Înțelegerea conceptului de integrare multi-omics (genomică, transcriptomică, epigenomică, proteomică).  
- Exersarea metodelor computaționale de integrare (corelații, PCA comună, network fusion).  
- Construirea unui mini-pipeline multi-omics pe date sintetice / reduse (SNP + expresie).  
- Identificarea markerilor comuni între straturi (cross-omics markers).  
- Înțelegerea aplicabilității în cancer, boli rare și farmacogenomică.

---

## Context
După ce în săptămânile 6–9 am studiat:

- rețele de gene (GCE)  
- vizualizare & diseasome  
- repurposing bazat pe rețele  
- modele ML (supervised / unsupervised)

acum vom integra **multiple tipuri de date biologice** într-un singur spațiu analitic.

Integrarea multi-omics este esențială pentru:
- **subtipare tumorală**,  
- **stratificarea pacienților**,  
- **descoperirea de biomarkeri**,  
- **identificarea țintelor terapeutice**.

Vom lucra pe un set simplificat (SNPs + expresie + opțional methylation) pentru a înțelege logica pipeline-ului.

---

## Hands-on

### 🔹 `demo01_multiomics_alignment.py`  
- încarcă 2 matrici (SNP × samples, Expression × samples),  
- normalizează fiecare layer,  
- aliniază probele comune,  
- construiește un feature matrix concat.pkl.

### 🔹 `demo02_pca_fusion.py`  
- aplică PCA pe date integrate,  
- compară PCA_snp, PCA_expr, PCA_joint,  
- vizualizează probele.

### 🔹 `demo03_correlation_integration.py`  
- calculează corelații SNP–gene,  
- construiește un grafic bipartit SNP–gene,  
- extrage perechi candidate pentru biomarkeri.

---

## Exerciții

### Exercițiul 1 — Feature Fusion & PCA  
- construiți matricea multi-omics (concat),  
- rulați PCA,  
- salvați încărcăturile și proiecțiile PCA.

### Exercițiul 2 — Cross-Omics Marker Discovery  
- identificați cele mai corelate perechi SNP–gene,  
- filtrați la threshold (|r| > 0.5),  
- exportați listele de candidați.

---

## Livrabile
1. `multiomics_concat_<handle>.csv`  
2. `pca_joint_<handle>.png`  
3. `snp_gene_pairs_<handle>.csv`  
4. `report_lab10_<handle>.pdf` (max 3 pagini)  
5. checklist PR completat

---

## Săptămâna următoare
- Deep learning pentru multi-omics (CNNs, autoencoders)  
- Graph embeddings pentru integrative drug discovery  
- Pregătirea proiectului final

---

## Competențe
- Preprocesarea și armonizarea datelor multi-omics.  
- Aplicarea PCA pe date integrate.  
- Descoperirea perechilor SNP–gene candidate la biomarkeri.  
- Înțelegerea pipeline-urilor multi-omics folosite în genomica clinică.

---

## Resurse
- TCGA-BRCA, TCGA-GBM – multi-omics atlas  
- iCluster, MOFA, OmicsIntegrator  
- HapMap & GTEx cross-omics correlation  
- Jolliffe – PCA  
- Broad Institute: Single-Cell + multi-omics pipelines


