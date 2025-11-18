# Săptămâna 6 — Gene Co-Expression Networks (GCEs)  
Autor: IrisDanila

## Date folosite
Am folosit un dataset sintetic de expresie genică creat pentru a simula date RNA-Seq din cancer mamar:
- **Fișier**: `data/work/IrisDanila/lab06/expression_matrix.csv`
- **Dimensiune inițială**: 220 gene × 20 eșantioane
- **Structură**: 4 module teoretice de câte 50 gene fiecare + 20 gene random background
- **Design**: Module cu patternuri distincte de co-expresie între eșantioane

Am rulat următoarele scripturi:
```bash
python create_expression_data.py
python demo01_corr_threshold.py
python ex01_gce_networks.py


Metrică de corelație
Metodă: Spearman correlation
Justificare:
Mai robustă la outlier-i decât Pearson
Nu presupune relații liniare între gene
Recomandată pentru date biologice cu distribuții complexe
Abs(correlation): Da - folosim valoarea absolută pentru a captura atât corelații pozitive cât și negative
Prag pentru adiacență
Threshold: 0.6
Interpretare: Două gene formează o muchie dacă |cor| ≥ 0.6
Tip matrice: Binară (0/1), neponderat
Preprocesare
Log-transformare: log₂(x + 1) - normalizează distribuția și stabilizează varianța
Filtrare varianță: Threshold = 0.5
Eliminează gene constitutive (expresie constantă)
Păstrează doar gene cu variație semnificativă între eșantioane


Reflecție: Cum diferă o rețea de co-expresie față de clustering-ul clasic?
Deși ambele metode grupează entități similare, există diferențe fundamentale în abordare și rezultate:

1. Tipul de relații captate
Clustering (Lab 5 - K-means, Hierarchical, DBSCAN):

Relații globale: Fiecare eșantion/genă e atribuit unui singur cluster
Partiționare: Datele sunt divizate în grupuri disjuncte (non-overlapping)
Exemplu Lab 5: Pacient aparține fie Cluster 0 (benign) fie Cluster 1 (malign)
Nu captează: Relații specifice pereche-pereche
Rețele de co-expresie (Lab 6):

Relații locale pereche-pereche: Fiecare pereche de gene are o muchie dacă corelația e ridicată
Overlap posibil: O genă poate fi în multiple module (în funcție de algoritm)
Exemplu Lab 6: Gene_M1_001 e conectată direct cu Gene_M1_002 (cor > 0.6)
Captează: Structura detaliată a interacțiunilor
2. Reprezentarea datelor
Clustering:

Input: Matrice de features (eșantioane × caracteristici)
Output: Vector de etichete (fiecare eșantion → cluster_id)
Vizualizare: Scatter plot cu culori per cluster (2D/3D după PCA)
Informație: Apartenență la grup, centroizi, silhouette score
Rețele:

Input: Matrice de corelație sau adiacență (gene × gene)
Output: Graf (noduri = gene, muchii = co-expresie)
Vizualizare: Rețea cu noduri și muchii, layout force-directed
Informație: Conectivitate, hub genes, module, centralitate


Concluzie
Clustering și rețelele de co-expresie abordează întrebări complementare:

Clustering: "Care entități sunt similare global?" → Grupuri omogene
Rețele: "Cine interacționează cu cine?" → Structură de relații
În Lab 5, am folosit clustering pentru a stratifica pacienți cu cancer mamar în subtipuri. În Lab 6, am construit rețele pentru a identifica module de gene co-exprimate.

Pentru bioinformatică integrativă:

Clustering pe eșantioane identifică ce diferențiază grupurile
Rețele de gene explică cum (mecanisme moleculare) apar diferențele