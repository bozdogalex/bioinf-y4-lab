
### Rezultatul rularii demo01_distance_matrix.py
python '/workspaces/bioinf-y4-lab/labs/04_phylogenetics/demo01_distance_matrix.py'
Sequences: ['NM_000546.6', 'NM_011640.3', 'NM_131327.2']
Distance matrix:
 [[0.         0.51194268 0.6691879 ]
 [0.51194268 0.         0.72599663]
 [0.6691879  0.72599663 0.        ]]






# Analiza Filogenetica - Genele p53 (TP53)


## Metode
* **Gene Utilizate:** Trei secvente de ARNm ale genei TP53/Trp53.
* **Algoritm Filogenetic:** Neighbor-Joining (NJ).
* **Distante:** Bazate pe alinierea secventelor (valori intre 0 si 1, indicand fractiunea de diferente).

## Rezultate si Matricea de Distante

Genele folosite au fost:
1. **NM_000546.6** — Homo sapiens TP53 mRNA
2. **NM_011640.3** — Mus musculus Trp53 mRNA
3. **NM_131327.2** — Danio rerio tp53 mRNA

**[INFO] Matricea de distante:**

| Taxon | SRR6024194.1 | SRR6024194.2 | SRR6024194.3 |
| :--- | :---: | :---: | :---: |
| **SRR6024194.1** | 0.000000 | | |
| **SRR6024194.2** | 0.786667 | 0.000000 | |
| **SRR6024194.3** | 0.746667 | 0.786667 | 0.000000 |

* **Distanta minima:** 0.746667 (intre SRR6024194.1 si SRR6024194.3).
* **Distanta maxima:** 0.786667 (intre 1 si 2, si 2 si 3).

## Arborele Filogenetic (ASCII Tree)
  _________________________________________________________ SRR6024194.1
 |
_|________________________________________________________________ SRR6024194.2
 |
 |_________________________________________________________ SRR6024194.3

 ## Interpretare Filogenetica

Arborele NJ indica faptul ca secventele **SRR6024194.1** si **SRR6024194.3** sunt perechea cu cea mai mica distanta genetica ($0.746667$), sugerand o divergenta mai recenta intre ele comparativ cu **SRR6024194.2**.

In contextul filogeniei speciilor (Mamifere grupate), rezultatul obtinut este neasteptat, deoarece ne-am fi asteptat ca ortologii de la Homo sapiens si Mus musculus sa fie cei mai apropiati. Aceasta poate indica:
1. O **etalonare neobisnuita** a secventelor (e.g., SRR6024194.1 = H. sapiens, SRR6024194.3 = M. musculus, SRR6024194.2 = D. rerio).
2. O **evolutie diferita a genei p53** (ex. rate de substitutie diferite sau constrangeri functionale) in regiunile analizate care nu reflecta fidel arborele speciilor.