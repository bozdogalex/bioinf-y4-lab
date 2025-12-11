### DEMO 01

Sequences: ['NM_000546.6', 'NM_011640.3', 'NM_131327.2']
Distance matrix:
 [[0.         0.51194268 0.6691879 ]
 [0.51194268 0.         0.72599663]
 [0.6691879  0.72599663 0.        ]]


### EX 01

[INFO] Toate secvențele au aceeași lungime -> folosesc p-distance (Hamming).
[OK] Arbore NJ salvat -> labs/04_phylogenetics/submissions/Razvann19/tree_Razvann19.nwk

=== NJ Tree (ASCII preview) ===
                                      , NG_017013.2_2
  ____________________________________|
 |                                    | NG_017013.2_1
 |
 |____________________________________ NM_000546.6_1
_|
 |                                   , NC_060941.1_2
 |___________________________________|
 |                                   | NC_060941.1_1
 |
 |_____________________________________________________________ NC_000017.11_1



### Am folosit trei secvențe TP53 (Homo sapiens), descărcate în laboratorul 1:

| Fișier | Descriere | Link |
|-------|-----------|------|
| `my_tp53.fa` | Secvență TP53 Homo sapiens descărcată prin NCBI Entrez | https://www.ncbi.nlm.nih.gov/nuccore |
| `nm000546.fa` | Transcript referință TP53 (NM_000546.6) — mRNA | https://www.ncbi.nlm.nih.gov/nuccore/NM_000546 |
| `my_tp53_small.fa` | Versiune scurtată a TP53 folosită pentru testare/aliniere rapidă | generată local din `my_tp53.fa` |

Fișierul multi-FASTA a fost salvat în `data/work/Razvann19/lab04/your_sequences.fasta`.

---

###  Reflectie

**Ce informații suplimentare oferă arborele filogenetic față de o simplă matrice de distanțe?**

O matrice de distanțe arată doar *cât de diferite* sunt secvențele între ele, dar nu indică *structura relațiilor lor evolutive*.  

Un arbore filogenetic oferă:

- Relațiile ierarhice dintre secvențe (cine este mai apropiat de cine)
- Un punct de vedere evolutiv asupra divergenței
- Topologia evoluției — ordinea ramificării
- O reprezentare vizuală intuitivă a clustering-ului secvențelor

Astfel, arborele nu doar cuantifică distanța genetică, ci **interpretează și organizează** aceste diferențe într-un model evolutiv coerent.

---


