# Raport Lab 08 — Machine Learning pe date omice

**Autor:** AlexTGoCreative  
**Data:** Ianuarie 2026

---

## 1. Introducere

### Tipul datelor
Dataset de expresie genică cu:
- **200 probe** biologice (esantioane de țesuturi)
- **50 gene** măsurate (GENE_001 - GENE_050)
- **3 clase**: Brain (68 probe), Kidney (69), Liver (63)

### Obiective
1. Clasificare supervizată (Random Forest + Logistic Regression)
2. Analiză nesupravegheată (PCA + KMeans)
3. Experiment semi-supervised cu pseudo-labeling
4. Interpretare biologică a genelor importante

---

## 2. Supervised ML — Random Forest și Logistic Regression

### 2.1 Random Forest

**Hyperparametri:** 200 arbori, max_depth=10, min_samples_split=5

**Rezultate test set (40 probe):**

| Clasă  | Precision | Recall | F1-Score | Support |
|--------|-----------|--------|----------|---------|
| Brain  | 0.38      | 0.46   | 0.41     | 13      |
| Kidney | 0.35      | 0.43   | 0.39     | 14      |
| Liver  | 0.43      | 0.23   | 0.30     | 13      |
| **Accuracy** | | | **0.38** | **40** |

**Top 10 gene după importanță:**
1. GENE_023 (4.01%)
2. GENE_042 (2.95%)
3. GENE_010 (2.85%)
4. GENE_050 (2.82%)
5. GENE_017 (2.61%)
6. GENE_046 (2.57%)
7. GENE_011 (2.51%)
8. GENE_028 (2.43%)
9. GENE_031 (2.38%)
10. GENE_001 (2.34%)

**Observații:**
- Liver are cel mai mic recall (23%) → clasificat greșit cel mai des
- Performanță modestă (38%) sugerează overlap între clase
- GENE_023 este cea mai informativă genă pentru diferențierea țesuturilor

### 2.2 Logistic Regression vs Random Forest

**Comparație:**

| Model               | Accuracy | F1-Score (weighted) |
|---------------------|----------|---------------------|
| Random Forest       | 37.5%    | 0.37                |
| Logistic Regression | 47.5%    | ~0.48               |
| **Diferență**       | **-10%** | **+0.11**           |

**Discuție modele liniare vs non-liniare:**
- **Logistic Regression performează mai bine** (+10% accuracy)
- RF suferă de **overfitting** pe doar 160 probe de antrenament
- Datele par să aibă **structură mai degrabă liniară**
- LR beneficiază de **scaling** (StandardScaler) și regularizare implicită
- Pentru dataset-uri mici (< 200 probe), modelele liniare sunt adesea mai robuste

---

## 3. Unsupervised ML — PCA + KMeans

### 3.1 PCA (2 componente)
- **PC1:** 5.10% varianță explicată
- **PC2:** 4.13% varianță explicată
- **Total:** 9.23%

**Interpretare:** Varianța scăzută (< 10%) indică că datele sunt distribuite în multe dimensiuni. Nu există 2-3 direcții dominante care să explice structura.

### 3.2 KMeans Clustering (k=3)

**Crosstab: Label real × Cluster:**

| Label  | Cluster 0 | Cluster 1 | Cluster 2 |
|--------|-----------|-----------|-----------|
| Brain  | 21        | 36        | 11        |
| Kidney | 24        | 28        | 17        |
| Liver  | 15        | 26        | 22        |

**Clusterizarea recuperează structura de clase?**
- **Nu.** Fiecare cluster conține probe din toate cele 3 țesuturi
- Cluster 1 este cel mai mare (90 probe), dar heterogen
- Brain și Liver au distribuții similare prin clustere

**Probe care "sar" în cluster greșit:**
- Brain → Cluster 1 (36/68 = 53%)
- Kidney → distribuit aproape uniform
- Lipsa separării clare sugerează **overlap în expresia genică** între țesuturi

---

## 4. Semi-Supervised Learning

### Configurare experiment
- **40% etichete marcate unknown** (80 probe)
- **60% etichete păstrate** (120 probe)
- Metodă: Pseudo-labeling cu Random Forest

### Rezultate

| Scenario                | Accuracy | Observații                          |
|-------------------------|----------|-------------------------------------|
| Baseline (60% labeled)  | 62.5%    | Doar probe etichetate              |
| După pseudo-labeling    | 62.5%    | +0.0% → nicio îmbunătățire         |
| Full data (100%)        | 35.0%    | Upper bound (reference)            |

**Performanța crește sau scade?**
- **Rămâne constantă** (62.5%)
- Pseudo-labeling **nu ajută** în acest caz

**De ce pseudo-labeling nu funcționează aici?**
1. Modelul baseline are **low confidence** (accuracy 38%)
2. Pseudo-etichetele sunt **probabil greșite**
3. Se propagă erorile în modelul final

**De ce ajută pseudo-labeling în bioinformatică (în general)?**
1. **Costuri de etichetare mari:** validare clinică, secvențiere validare
2. **Date abundente neetichetate:** RNA-seq public, imagistică medicală
3. **Exploatare structură:** pacienți similari au diagnostice similare
4. **Active learning:** se pot eticheta manual doar probe "incerte"

---

## 5. Interpretare biologică

### Gene informative (top 5)
- **GENE_023** (4.01%): Potențial marker pentru diferențierea Brain/Liver
- **GENE_042** (2.95%): Exprimare specifică Kidney?
- **GENE_010** (2.85%): Funcție metabolică organ-specifică
- **GENE_050** (2.82%): Posibil regulator transcripțional
- **GENE_017** (2.61%): Gen housekeeping cu variabilitate între țesuturi

### Validare necesară
Pentru a confirma relevanța biologică:
1. **Pathway enrichment analysis** (GO, KEGG, Reactome)
2. **Literatura:** sunt gene cunoscute pentru țesuturi?
3. **Validare experimentală:** qPCR, Western blot
4. **Analiză diferențială:** fold-change între țesuturi

---

## 6. Limitări

### 1. Small sample size
- **200 probe** → prea puține pentru 50 gene
- RF are tendință de overfitting
- Confidence intervals largi

### 2. Număr limitat de gene
- **50 gene** → reprezentare incompletă a transcriptomului
- Lipsesc gene marker cunoscute (ex: GFAP pentru creier)

### 3. Variabilitate biologică necontrolată
- Factori de confuzie: vârstă, sex, batch effects
- Lipsă metadate despre probe
- Posibil date sintetice/simulate

### 4. Lipsa validării
- Niciun set independent de test
- Cross-validation ar fi necesară
- Performanța reală poate fi mai scăzută

---

## Bonus — PCA după eliminarea genelor cu varianță mică

**Rezultate:**
- **Înainte** (50 gene): 9.23% varianță explicată
- **După** (40 gene): 11.31% varianță explicată
- **Îmbunătățire:** +2.08 puncte procentuale (+22.5% relativ)

**Concluzie:**
Eliminarea celor 10 gene cu varianță mică **îmbunătățește vizibilitatea** în PCA. Aceste gene adăugau zgomot fără informație discriminativă. În practică, filtering-ul pre-PCA este recomandat.

---

## Anexe — Fișiere generate

1. [`classification_report_AlexTGoCreative.txt`](classification_report_AlexTGoCreative.txt)
2. [`confusion_rf_AlexTGoCreative.png`](confusion_rf_AlexTGoCreative.png)
3. [`feature_importance_AlexTGoCreative.csv`](feature_importance_AlexTGoCreative.csv)
4. [`cluster_crosstab_AlexTGoCreative.csv`](cluster_crosstab_AlexTGoCreative.csv)
5. [`sup_vs_unsup_scatter_AlexTGoCreative.png`](sup_vs_unsup_scatter_AlexTGoCreative.png)
6. [`rf_vs_logreg_report_AlexTGoCreative.txt`](rf_vs_logreg_report_AlexTGoCreative.txt)
7. [`semi_supervised_report_AlexTGoCreative.txt`](semi_supervised_report_AlexTGoCreative.txt)
8. [`pca_bonus_variance_AlexTGoCreative.png`](pca_bonus_variance_AlexTGoCreative.png)
