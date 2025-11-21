 # Lab 04 — Phylogenetics — Notes (olezhka0809)

## Results of running demo02_distance_matrix.py
Sequences: ['NM_000546.6', 'NM_011640.3', 'NM_131327.2']
Distance matrix:
 [[0.         0.51194268 0.6691879 ]
 [0.51194268 0.         0.72599663]
 [0.6691879  0.72599663 0.        ]]



## Secvențele FASTA folosite

Pentru construirea arborelui Neighbor-Joining am folosit un fișier multi-FASTA cu 5 secvențe proteice corespunzătoare proteinei p53 umane (`TP53`), forma canonică și patru izoforme:

- `sp|P04637.4|P53_HUMAN` — Cellular tumor antigen p53 (forma canonică, UniProt)
- `NP_001263628.1` — p53 isoform l [Homo sapiens]
- `NP_001263627.1` — p53 isoform k [Homo sapiens]
- `NP_001263626.1` — p53 isoform j [Homo sapiens]
- `NP_001263625.1` — p53 isoform i [Homo sapiens]

Fișierul multi-FASTA este salvat în:

`/workspaces/bioinf-y4-lab/data/work/olezhka0809/lab04/your_sequences.fasta`

Acest fișier a fost creat pornind de la secvența TP53 descărcată și folosită în laboratorul 02, la care am adăugat izoformele suplimentare ale proteinei p53 din baza de date NCBI/UniProt.

## Reflecție: arbore filogenetic vs. matrice de distanțe

O simplă matrice de distanțe arată doar „cât de diferite” sunt secvențele două câte două, dar nu spune nimic clar despre structura de grupare sau despre relațiile ierarhice dintre ele. Arborele filogenetic transformă aceste distanțe într-o reprezentare structurală, unde se văd direct clusterele de secvențe care sunt mai apropiate între ele și modul în care aceste grupuri se despart de la un „strămoș” comun. În plus, lungimea ramurilor din arbore oferă o intuiție vizuală despre cât de mari sunt diferențele, nu doar numeric, ci și în contextul celorlalte secvențe. Practic, arborele condensează informația din matrice într-o formă care este mult mai ușor de interpretat biologic (de exemplu: care izoforme sunt cele mai apropiate de forma canonică p53 și care sunt mai divergente).
