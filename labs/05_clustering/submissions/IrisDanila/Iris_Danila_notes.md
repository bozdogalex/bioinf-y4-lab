# Săptămâna 5 — Clustering în Bioinformatică  
Autor: IrisDanila

## Date folosite
Am folosit dataset-ul **Breast Cancer Wisconsin Diagnostic (WDBC)** de pe UCI Repository:
- **Sursă**: https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data
- **Dimensiune**: 569 eșantioane, 30 features (caracteristici nuclei celulari)
- **Clase**: 
  - Malignant (M): 212 cazuri (37.3%)
  - Benign (B): 357 cazuri (62.7%)
- **Features**: Caracteristici calculate din imagini digitale ale aspiratelor cu ac fin (FNA) ale maselor mamare

Am rulat următoarele scripturi:
```bash
python demo01_k_means.py
python ex01_clustering.py


Rezultate și observații
Preprocesare și Standardizare
Dataset standardizat: 569 eșantioane × 30 features
Toate features au fost scalate la mean=0, std=1
Esențial pentru metodele bazate pe distanță (K-means, DBSCAN)
1. Hierarchical Clustering (Average Linkage)
Rezultate:

Dendrogramă generată cu succes
Metoda: Average linkage
Vizualizare truncată la ultimele 30 clustere pentru claritate
Observații:

Dendrograma arată o structură ierarhică clară
Se pot identifica potențial 2 clustere majore (corespunzând malign/benign)
Lungimea ramurilor indică distanța între clustere
Utilă pentru a determina numărul optim de clustere
2. K-means Clustering (K=2)
Rezultate:

Code
Cluster distribution:
  Cluster 0: 375 eșantioane
  Cluster 1: 194 eșantioane

PCA Explained Variance:
  PC1: 44.27%
  PC2: 18.97%
  Total: 63.24%
Observații:

K-means a identificat 2 clustere distincte
Cluster 1 (194 eșantioane) se potrivește aproximativ cu cazurile maligne (212)
Cluster 0 (375 eșantioane) conține cazuri benigne + unele maligne
PCA păstrează ~63% din variația originală în 2D
Separarea vizuală în plot este destul de clară
Potrivire bună între clustere și diagnostic real
Performanță:

Relativ bună pentru un algoritm nesupervizat
Asumă clustere sferice și de dimensiune similară
Convergență rapidă și reproductibilă (random_state=42)
3. DBSCAN Clustering
Rezultate:

Number of clusters: 1
Number of noise points: 550
Cluster distribution:
  Cluster 0: 19 eșantioane
  Noise (-1): 550 eșantioane
Observații:

DBSCAN cu parametrii eps=1.5, min_samples=5 a eșuat să găsească structură semnificativă
96.7% din date clasificate ca "noise"
Doar 19 eșantioane într-un singur cluster dens
Concluzie: Parametrii nu sunt optimali pentru acest dataset

Reflecție: Cum se compară clustering-ul cu arborii filogenetici în descoperirea relațiilor biologice?
Deși ambele metode grupează entități similare, există diferențe fundamentale în aplicare și interpretare:

1. Scopuri diferite
Arbori filogenetici:

Scop: Reconstrucția istoriei evolutive și relațiilor de rudenie
Focus: Relații ancestrale și ordinea ramificărilor în timp
Direcție: Arată cine provine de la cine (rădăcină → frunze)
Exemplu Lab 4: Variantele TP53 grupate după similaritatea secvențelor, sugerând splicing alternativ comun
Clustering:

Scop: Descoperirea grupurilor naturale în date, fără presupuneri despre istorie
Focus: Similaritate/proximitate în spațiul caracteristicilor
Direcție: Nu există direcție temporală sau cauzală
Exemplu Lab 5: Pacienți grupați după caracteristicile tumorilor, pentru diagnostic



Concluzie
Arborii filogenetici și clustering-ul sunt două fațete ale aceluiași obiectiv - descoperirea relațiilor biologice - dar din perspective diferite:

Filogeniile presupun un model evolutiv și caută ISTORIA relațiilor
Clustering-ul explorează date fără presupuneri și caută GRUPĂRI actuale
În Lab 4, am folosit phylogenetics pentru a înțelege relațiile evolutive ale variantelor TP53. În Lab 5, am folosit clustering pentru a stratifica pacienți cu cancer mamar.