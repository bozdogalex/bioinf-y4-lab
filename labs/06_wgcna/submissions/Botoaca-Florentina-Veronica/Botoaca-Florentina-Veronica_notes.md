# Note pentru Laboratorul 6 - Gene Co-Expression Networks

## Parametrii utilizați
- **Metrică de corelație**: Spearman (mai robust la outliers față de Pearson)
- **Prag de corelație**: 0.6 (|cor| ≥ 0.6 pentru a crea o muchie)
- **Prag de varianță**: 0.5 (pentru filtrarea genelor cu expresie constantă)
- **Folosim valori absolute**: Da (corelațiile negative puternice sunt la fel de relevante)
- **Graf neorientat**: Da
- **Algoritm de detectare a modulelor**: Louvain (cu fallback la Greedy Modularity)

## Observații practice
- Am folosit log2(x+1) pentru a normaliza distribuția datelor RNA-Seq
- Filtrarea genelor cu varianță mică este crucială pentru a reduce zgomotul
- Pragul de 0.6 a fost ales pentru a echilibra între specificitate și sensibilitate
- Louvain algoritmul este preferat pentru detectarea modulelor datorită eficienței



============================================================
LABORATOR 6: GENE CO-EXPRESSION NETWORKS
Handle: Botoaca-Florentina-Veronica
============================================================
Matrice încărcată: 200 gene, 50 probe
Valori min/max: 4.74 / 15.29

=== PREPROCESARE DATES ===
Date intrare: 200 gene
✓ Transformare log2(x+1) aplicată

Statistici varianță:
  Min:    0.00241835
  Medie:  0.03950796
  Max:    0.07790595
  Prag (50%): 0.04656002

Distribuție varianțe:
    0%: 0.00241835
   25%: 0.02214304
   50%: 0.04656002
   75%: 0.05895955
   95%: 0.06939434
  100%: 0.07790595

✓ Gene care trec pragul: 100 / 200
✓ Gene rămase după filtrare: 100

=== CALCUL CORELAȚII ===
Metodă: spearman
✓ Folosesc valorile absolute ale corelației
✓ Matrice corelație: 100×100
  Valori: [0.000, 1.000]

=== CONSTRUIRE ADIACENȚĂ ===
Prag pentru muchie: 0.6
✓ Matrice binară
✓ Muchii: 1650
✓ Densitate graf: 0.333333

=== CONSTRUIRE GRAF ===
✓ Graf neorientat creat
✓ Noduri finale: 100
✓ Muchii finale: 1650

=== DETECTARE MODULE ===
✓ Algoritm: Louvain
✓ Module detectate: 3
✓ Dimensiuni module:
  Min: 30 gene
  Max: 40 gene
  Medie: 33.3 gene
✓ Fișier salvat: /workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/Botoaca-Florentina-Veronica/modules_Botoaca-Florentina-Veronica.csv
✓ Total gene în module: 100

============================================================
REZUMAT FINAL:
- Gene inițiale: 200
- Gene după filtrare: 100
- Noduri în graf: 100
- Muchii în graf: 1650
- Module detectate: 3
============================================================