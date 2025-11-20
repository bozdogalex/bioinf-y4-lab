# Demo 01 - Correlation Threshold
```
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
```

# Exercitiu 1 - GCE Networks
```
Grafic creat cu 83 noduri și 974 muchii.
S-au detectat 8 module.
Am salvat mapping-ul gene→modul în: labs/06_wgcna/submissions/numipasaaa/modules_numipasaaa.csv
```

# Livrabile

- ce metrică de corelație și ce prag ați folosit,  
       - **Metrica de corelatie**: Spearman
       - **Pragul de adiacenta**: 0.6
       - **Pragul de variatie**: 0.5
- o scurtă reflecție: **Cum diferă o rețea de co-expresie față de clustering-ul clasic (Lab 5)?**
       - **Clustering-ul clasic**: Obiectivul este **partitionarea**. Datele sunt fortate intr-un numar specific de grupuri (clustere) discrete, unde fiecare gena apartine unui singur cluster. Rezultatul este o **lista de membri** pentru fiecare grup.
       - **Reteaua de co-expresie**: Obiectivul este **definirea conectivitatii**. Se calculeaza o matrice de adiacenta (conexiuni) intre toate perechile de gene. Desi genele sunt grupate in module, rezultatul principal este **arhitectura retelei**.
       - O retea ofera informatii structurale pe care clustering-ul clasic **nu** le poate furniza:
              - **Conectivitate intra-modulara (Hub Genes)**: Clustering-ul trateaza toate genele dintr-un grup ca fiind egale. Reteaua identifica **genele "hub"**, nodurile cele mai inalt conectate din cadrul unui modul.
              - **Relatii inter-modulare**: Clustering-ul produce grupuri izolate. O retea poate arata **conexiunile dintre module**, permitand modelarea interactiunilor la un nivel superior (retele de module).
              - **Masuratori topologice**: Reteaua permite calcularea unor metrici specifice (ex. centralitate, coeficient de clustering) care descriu importanta relativa a fiecarei gene in intregul sistem, nu doar apartenenta sa la un grup.
       - Clustering-ul doar **grupeaza** genele. O retea de co-expresie **modeleaza relatiile** dintre ele, identificand atat grupurile (modulele), cat si arhitectura de reglare din cadrul si dintre acestea.
