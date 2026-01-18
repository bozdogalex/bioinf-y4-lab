# Assignment 5 - Clustering și Filogenetică
## Analiza Comparativă între Clustering și Arbore Filogenetic

**Student:** AlexTGoCreative  
**Data:** Decembrie 2024  
**Dataset:** 10 secvențe proteice p53 din multiple specii

---

## 1. Rezumat Metodologie

### 1.1 Dataset și Preprocesare
Am reutilizat dataset-ul multi-FASTA din Lab 4, conținând 10 secvențe proteice ale genei p53 de la:
- Homo sapiens (ASE05898.1)
- Mus musculus (NP_001120705.1)
- Rattus norvegicus (sp|P10361.1|P53_RAT)
- Danio rerio (NP_001315517.1)
- Xenopus tropicalis (NP_001001903.1)
- Gallus gallus (NP_990595.1)
- Bos taurus (NP_776626.1)
- Sus scrofa (NP_001376147.1, NP_998989.3)
- Equus caballus (XP_077825621.1)

**Reprezentare:** K-mer frequency vectors (k=3) → matrice 10×1261 caracteristici normalizate

### 1.2 Metode de Clustering Aplicate

#### **A. K-means Clustering**
- Testat cu k = 2, 3, 4 clustere
- **Rezultate (k=2 - cel mai bun):**
  - Silhouette Score: **0.2359** (separare moderată)
  - Davies-Bouldin Index: 0.5631 (compactitate bună)
  - Calinski-Harabasz Index: 2.4261

**Cluster 0:** 9 secvențe (toate mamifere + păsări + amfibian + pește)  
**Cluster 1:** 1 secvență (ASE05898.1 - Homo sapiens partial)

#### **B. Hierarchical Clustering (Average Linkage)**
- Testat cu k = 2, 3, 4 clustere
- **Rezultate (k=2 - cel mai bun):**
  - Silhouette Score: **0.2359** (identic cu K-means)
  - Rezultat aproape identic cu K-means pentru k=2

#### **C. DBSCAN**
- Testat cu eps = 0.3, 0.5, 0.7 și min_samples = 2
- **Rezultate:** Un singur cluster pentru toate parametrii
- **Interpretare:** Dataset-ul este prea omogen pentru separare bazată pe densitate

---

## 2. Comparație cu Arborele Filogenetic

### 2.1 Observații Principale

**Concordanță parțială:** Clustering-ul a separat doar o singură secvență (ASE05898.1 - Homo sapiens partial) de restul grupului.

**Discrepanță majoră:** Arborele filogenetic arată multiple clade-uri clare:
- Clade mamifere (șoareci, șobolani, bovine, porcine)
- Clade păsări/reptile (Gallus)
- Clade amfibian (Xenopus)
- Clade pești (Danio rerio)

### 2.2 Analiza Discrepanțelor

| Aspect | Clustering | Filogenetică |
|--------|-----------|--------------|
| **Bază** | Similaritate k-mer (secvență locală) | Distanță evolutivă (istorie comună) |
| **Rezultat** | 2 grupuri (1 secvență vs. 9 secvențe) | Multiple clade-uri evolutive |
| **Interpretare** | Omogenitate funcțională | Diversitate evolutivă |

### 2.3 Ipoteze Biologice

**1. Conservare funcțională extremă a p53:**
- Proteina p53 este un "gardian al genomului" cu funcție critică
- Presiune selecțivă puternică → secvențe foarte conservate
- K-mer similarities ridicate chiar între specii distante evolutiv

**2. Secvența ASE05898.1 (outlier):**
- Marcată ca "partial" → posibil fragmentată sau incompletă
- Lipsa unor domenii conservate → separare în cluster distinct
- Nu reflectă divergență biologică reală, ci artefact de date

**3. Evoluție convergentă vs. Omologie:**
- Clustering capturează similarități funcționale (domenii active)
- Filogenetica capturează istorie evolutivă (mutații neutre + active)

---

## 3. Evaluare Metode de Clustering

### 3.1 Performance Metrics

| Metodă | Silhouette | Davies-Bouldin | Calinski-Harabasz | Interpretare |
|--------|-----------|----------------|-------------------|--------------|
| **K-means (k=2)** | 0.2359 | 0.5631 | 2.4261 | Cel mai bun |
| **Hierarchical (k=2)** | 0.2359 | 0.5631 | 2.4261 | Identic K-means |
| **K-means (k=3)** | 0.1520 | 1.1619 | 2.4553 | Separare slabă |
| **DBSCAN (toate)** | 0.0000 | N/A | N/A | Nerelevant |

**Concluzie:** K-means și Hierarchical cu k=2 oferă cel mai bun compromis, dar Silhouette Score scăzut (0.24) indică structură de clustering slabă în date.

### 3.2 Diferențe între Metode

- **K-means vs. Hierarchical:** Rezultate identice pentru k=2 (dataset mic, omogen)
- **DBSCAN:** Eșuează să detecteze structură → toate secvențele în același cluster
- **Interpretare:** Dataset-ul nu are separare naturală în "clustere" bazate pe densitate

---

## 4. Interpretări Biologice și Funcționale

### 4.1 Perspective Funcționale

**Conservare domenii funcționale:**
- Domeniul de legare ADN (DNA-binding domain)
- Domeniul de tetramerizare
- Domeniul de transactivare

Toate acestea sunt extrem de conservate → similaritate mare în k-mer space

### 4.2 Aplicații Practice

**1. Identificarea subtipurilor de boală:**
- Mutații p53 în cancer → necesită analiză mai fină decât k-mer clustering
- Combinarea cu date structurale (3D) ar putea îmbunătăți clustering-ul

**2. Gene families și evoluție funcțională:**
- Pentru gene families mai diverse, clustering-ul poate identifica subfamilii funcționale
- P53 este prea conservat pentru această aplicație

**3. Drug design și repurposing:**
- Clustering poate identifica regiuni conservate pentru țintire terapeutică
- Filogenetica ajută la predicția efectelor cross-species

---

## 5. Concluzii

### 5.1 Răspunsuri la Întrebările de Reflecție

**Q1: Cum se aliniază rezultatele clustering-ului cu arborele filogenetic?**
- **Răspuns:** Aliniere slabă. Clustering-ul produce o separare minimă (9 vs. 1), în timp ce filogenetica arată multiple clade-uri distincte. Acest dezacord reflectă diferența între similaritate funcțională (clustering) și istorie evolutivă (filogenetică).

**Q2: Ce explică diferențele observate?**
- **Răspuns:** 
  1. **Conservare funcțională extremă:** p53 este vital → presiune selecțivă puternică
  2. **Reprezentare k-mer:** Capturează similarități locale, nu mutații sinonime
  3. **Artefact de date:** ASE05898.1 partial → outlier artificial

**Q3: Cum poate fi utilă combinarea clustering-ului cu filogenetica?**
- **Răspuns:** 
  - **Complementaritate:** Clustering identifică grupuri funcționale, filogenetica istorie evolutivă
  - **Validare:** Discrepanțe pot semnala evenimente evolutive interesante (gene duplication, horizontal transfer)
  - **Drug design:** Clustering pentru ținte conservate + filogenetică pentru specificitate

### 5.2 Limitări și Îmbunătățiri

**Limitări:**
- Dataset mic (n=10) → putere statistică redusă
- K-mer representation simplă → nu capturează structură 3D
- Secvență partial (ASE05898.1) → zgomot în date

**Îmbunătățiri posibile:**
- Alignment-based features (conserved domains)
- Structural features (predicted 3D structure)
- Functional annotations (GO terms, pathways)
- Larger dataset (mai multe specii + isoforms)

---

## 6. Livrabile Generate

### 6.1 Cod Python
- `clustering_phylogeny_analysis.py` - Script complet cu 3 metode de clustering

### 6.2 Vizualizări
- `dendrogram.png` - Dendrogramă hierarchical clustering
- `pca_kmeans.png` - Proiecție PCA cu clustere K-means
- `pca_hierarchical.png` - Proiecție PCA cu clustere hierarchical
- `phylo_tree_annotated.png` - Arbore filogenetic adnotat cu clustere

### 6.3 Date
- `cluster_assignments.csv` - Asignări clustere pentru fiecare secvență
- `analysis_report.txt` - Raport detaliat metrics

---

**Notă finală:** Acest exercițiu demonstrează că clustering-ul și filogenetica oferă perspective complementare, nu identice. Pentru gene foarte conservate ca p53, filogenetica este mai informativă decât clustering-ul bazat pe k-mers. Însă combinarea ambelor abordări oferă o înțelegere mai completă a relațiilor evolutive și funcționale.
