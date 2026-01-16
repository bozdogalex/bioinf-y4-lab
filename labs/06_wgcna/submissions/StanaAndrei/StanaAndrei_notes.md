demo01:
Matrice expresie (toy):
       Sample1  Sample2  Sample3  Sample4  Sample5
GeneA        5        4        6        5        4
GeneB        3        3        2        4        3
GeneC        8        9        7       10        8

Matrice corelație (Spearman):
          GeneA     GeneB     GeneC
GeneA  1.000000 -0.353553 -0.432590
GeneB -0.353553  1.000000  0.917663
GeneC -0.432590  0.917663  1.000000

Matrice adiacență cu prag 0.7:
       GeneA  GeneB  GeneC
GeneA      0      0      0
GeneB      0      0      1
GeneC      0      1      0

---------------------------------
ex01:

Matrice inițială: 20007 gene, 8444 probe.

---------------------------------

Clustering-ul clasic (ex: K-means, ierarhic) se concentrează pe a pune genele în grupuri , în timp ce rețeaua de co-expresie se concentrează pe a modela relațiile dintre perechile de gene folosind grafuri.