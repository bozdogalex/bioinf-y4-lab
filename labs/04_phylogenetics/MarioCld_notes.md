#rulare demo 1
root@codespaces-f90576:/workspaces/bioinf-y4-lab# /usr/local/bin/python /workspaces/bioinf-y4-lab/labs/04_phylogenetics/demo01_distance_matrix.py
Sequences: ['NM_000546.6', 'NM_011640.3', 'NM_131327.2']
Distance matrix:
 [[0.         0.51194268 0.6691879 ]
 [0.51194268 0.         0.72599663]
 [0.6691879  0.72599663 0.        ]]

#rulare exercitiu 5
root@codespaces-f90576:/workspaces/bioinf-y4-lab# /usr/local/bin/python /workspaces/bioinf-y4-lab/labs/04_phylogenetics/submissions/MarioCld/ex05_phylo_tree.py

Matricea de distanțe:
seq1_Ecoli  0.000000
seq2_Salmonella 0.037037    0.000000
seq3_Shigella   0.037037    0.055556    0.000000
    seq1_Ecoli  seq2_Salmonella seq3_Shigella

Arborele a fost salvat în: labs/04_phylogenetics/submissions/MarioCld/tree_MarioCld.nwk

## Secvențele FASTA folosite
Am folosit 3 secvențe ADN sintetice generate pentru testare:
>seq1_Ecoli
ATGCGTACGTAGCTAGCTAGCTAGCTAAGCTAGCTGATCGTAGCTAGCTGATCG

>seq2_Salmonella
ATGCGTACGTAGCTAGCTGGCTAGCTAAGCTAGCTGATCGTAGCTAGCTGATCA

>seq3_Shigella
ATGCGTACGTAGATAGCTAGCTAGCTAAGCTAGCTGATCGTAGCTAGCTGATCT

Matricea de distanțe arată doar cât de diferite sunt perechile de secvențe între ele.  
Arborele filogenetic, oferă informații ierarhice despre relațiile evolutive (cum sunt grupate secvențele și care este strămoșul comun cel mai apropiat)