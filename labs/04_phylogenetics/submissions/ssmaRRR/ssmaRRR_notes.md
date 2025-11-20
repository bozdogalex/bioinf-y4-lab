am rulat demo:
    Sequences: ['NM_000546.6', 'NM_011640.3', 'NM_131327.2']
    Distance matrix:
    [[0.         0.51194268 0.6691879 ]
    [0.51194268 0.         0.72599663]
    [0.6691879  0.72599663 0.        ]]

Am rulat ex1:
     , seq1
    |
    |________________________________________________________________________ seq2
    _|
    |________________________________________________________________________ seq3
    |
    |________________________________________________________________________ seq4

Am folosit aceasta secventa fasta:
    >seq1
    ATGCTAGCTAGCTACGATCGATCGATCGATCGATCG
    >seq2
    ATGCTAGCTAGCTACGATGGATCGATCGATCGATCG
    >seq3
    ATGCTAGCTAGCTACGATCGATCGATTGATCGATCG
    >seq4
    ATGCTAGCTAGCTACGATCGATCGATCGATAGATCG

Arborele filogenetic transformă matricea de distanțe într-o structură ierarhică, 
în care se văd clar grupele de secvențe și relațiile de „rudenie” (clade, secvențe surori). 
În loc de simple valori numerice pereche, arborele oferă o ipoteză despre istoricul comun și 
ordinea divergențelor dintre secvențe. Lungimile ramurilor indică vizual cantitatea de schimbări evolutive, 
făcând mai ușoară identificarea clusterelor și a secvențelor atipice.