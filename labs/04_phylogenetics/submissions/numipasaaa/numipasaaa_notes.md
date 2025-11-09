# Demo 1 - Distance Matrix
```
Sequences: ['NM_000546.6', 'NM_011640.3', 'NM_131327.2']
Distance matrix:
 [[0.         0.51194268 0.6691879 ]
 [0.51194268 0.         0.72599663]
 [0.6691879  0.72599663 0.        ]]
```

# Exercitiu 1 - Phylo NJ
- ce secvențe FASTA ați folosit (link/descriere):
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000932.1,NC_001320.1,NC_001666.2,NC_007898.3,NC_001879.2&rettype=fasta&retmode=text"

- o scurtă reflecție: **Ce informații suplimentare oferă arborele filogenetic față de o simplă matrice de distanțe?**
    - O matrice de distante ofera date numerice brute, cuantificand divergenta dintre perechi de secvente, fara a explica structura acestor diferente.
    - Un arbore filogenetic adauga mai multe informatii suplimentare esentiale, precum:
        - Topologia (relatiile de inrudire): Grupeaza secventele intr-o ierarhie. Arata care secvente formeaza clade si care sunt mai apropiate evolutiv una de alta.
        - Stramosii comuni (nodurile): Identifica punctele de divergenta ipotetice din trecut.
        - Cronologia relativa: Prin lungimea ramurilor si pozitia nodurilor (de la radacina spre varfuri), arborele ilustreaza ordinea temporala a evenimentelor de speciatie sau duplicare genica.
    - Pe scurt, arborele filogenetic transforma datele cantitative (distantele) intr-un model istoric calitativ (evolutia).

- output:

```
  _________________________________________________________ NC_000932.1
 |
 |          ____________________________________________________ NC_001666.2
 |_________|
_|         |______________________________________________________ NC_001320.1
 |
 |_________________________________________________________ NC_007898.3
 |
 |_________________________________________________________ NC_001879.2

```

