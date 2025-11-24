# Lab 06 — Gene Co-Expression Networks (GCE)
## Date și preprocesare

- Date: GSE115469 (single-cell RNA-Seq, PBMC)
- Fișier brut: `GSE115469_Data.csv`
- Am construit o matrice de expresie cu:
  - rânduri = gene
  - coloane = probe / celule
- Pentru a reduce dimensiunea și a nu depăși resursele din Codespaces:
  - am selectat **1000 de probe** (coloane)
  - am păstrat **5000 de gene** (ex. primele sau cele mai expresate)

Pași de preprocesare:
- Transformare **log2(x + 1)** pe toate valorile de expresie.
- Calcul varianță pe fiecare genă și filtrare cu:
  - `VARIANCE_THRESHOLD = 0.5`
  - după filtrare au rămas **33 de gene** din 5000 (cele mai variabile).

## Corelație și adiacență

- Metrică de corelație: **Spearman**
  - motiv: robustă la distribuții non-normale și la outlieri, potrivită pentru date de expresie RNA-Seq.
- Am folosit **valoarea absolută a corelației**:
  - `USE_ABS_CORR = True` → atât corelațiile pozitive, cât și cele negative puternice sunt considerate legături.
- Matrice de corelație obținută: **33 × 33**.

Construirea adiacenței:
- Prag: `ADJ_THRESHOLD = 0.6`
- `weighted = True` → muchiile păstrează valoarea corelației (nu doar 0/1).
- Inițial: 33 noduri și ~14 muchii.
- Am eliminat nodurile izolate (fără muchii), rezultând:
  - **14 noduri**, **14 muchii** în graful final.

## Rețeaua de co-expresie și modulele

- Graf construit cu `networkx`, neorientat:
  - `MAKE_UNDIRECTED = True`
- Algoritm de detectare module:
  - **Louvain** (`community-louvain`) — maximizare modularity.
- Rezultate:
  - **4 module** de co-expresie detectate.
  - mapping gene → modul salvat în:
    - `labs/06_wgcna/submissions/olezhka0809/modules_olezhka0809.csv`

## Observații

- Filtrarea puternică pe varianță (threshold = 0.5) reduce foarte mult numărul de gene (33/5000), dar păstrează doar genele cu profil de expresie foarte variabil, deci mai informative pentru co-expresie.
- Pragul de corelație destul de strict (`0.6`) produce o rețea relativ rară (14 muchii), dar modulele detectate sunt bine separate.
- Faptul că Louvain a detectat **4 module** sugerează existența unor grupuri de gene cu expresie co-variată, posibil legate de:
  - subtipuri celulare diferite,
  - căi biologice distincte,
  - sau stări funcționale separate (activare/imunitate etc.).
