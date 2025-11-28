# Raport Analiză Clustering - Cancer Mamar
Această analiză aplică metode de clustering pe dataset-ul WDBC (Wisconsin Diagnostic Breast Cancer) pentru a identifica pattern-uri în datele de cancer mamar.

### Performanța Metodelor
| Metodă | Adjusted Rand Score | Număr Clustere | Observații |
|--------|---------------------|----------------|------------|
| **K-means** | **0.654** | 2 | Performanță bună, corelează bine cu diagnosticul real |
| **DBSCAN** | -0.023 | 1 cluster + 550 outliers | Performanță slabă, prea mulți outliers |

### Distribuția Datelor
- **Total probe**: 569
- **Diagnostic Benign (B)**: 357 (62.7%)
- **Diagnostic Malign (M)**: 212 (37.3%)

## Analiza Comparativă

###  K-means (Recomandat)
- **ARI: 0.654** - corelare semnificativă cu diagnosticul real
- Identifică corect 2 clustere principale
- Potrivit pentru acest tip de date medicale

###  DBSCAN (Nerecomandat)
- **ARI: -0.023** - fără corelare cu diagnosticul real
- **550 outliers** (96.6% din date!) - parametri nepotriviți
- Nu reușește să identifice structuri clare


### Clusterele K-means corespund probabil:
- **Cluster 0**: Probe cu caracteristici de țesut benign
- **Cluster 1**: Probe cu caracteristici de țesut malign

### Limitări observate:
- Suprapunere între clustere (ARI < 1.0)
- Nevoie de feature engineering suplimentar
- DBSCAN nepotrivit pentru acest tip de date

### Metoda optimă pentru acest dataset:
**K-means cu K=2** oferă cele mai bune rezultate și este ușor de interpretat clinic.

##  Concluzii

Clustering-ul K-means demonstrează utilitatea în analiza datelor medicale de cancer, oferind o separare semnificativă între probele benigne și maligne. DBSCAN necesită tuning extensiv pentru a fi util pe acest tip de date.

**Recomandare finală**: K-means este metoda cea mai potrivită pentru analiza clustering pe dataset-ul de cancer mamar.

---
