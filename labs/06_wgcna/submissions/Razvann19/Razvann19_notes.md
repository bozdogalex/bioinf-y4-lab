demo 1:
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


## Parametri folosiți

- Metrica de corelație: **Pearson** (mai rapidă și potrivită pentru subsetul de gene standardizat, după log2(x+1)).
- Log-transformare: **log2(x + 1)** aplicată pe matricea de expresie RNA-Seq.
- Filtrare: am păstrat **cele mai variabile 4000 de gene** (din primele **8000 de rânduri** ale fișierului).
- Prag pentru adiacență: **|cor| ≥ 0.6**, pentru a păstra doar conexiunile puternice.
- Rețea neorientată; module detectate cu algoritmul **Louvain**, cu fallback pe greedy modularity.

## Notă privind fișierul de expresie
Fișierul original GSE115469 este foarte mare (483 MB), iar calculul integral al matricei de corelație pe toate genele depășește memoria disponibilă în mediul de lucru, ducând la întreruperea procesului (“terminated”).  
De aceea, am folosit primele **8000 de rânduri** din CSV și am selectat apoi **top 4000 gene** după varianță.  


## Reflecție — rețea de co-expresie vs. clustering clasic (Lab 5)

În Lab 5, clustering-ul (Hierarchical, K-means, DBSCAN) împarte genele sau probele în **clustere globale**, fiecare element aparținând unui singur cluster.  
În rețelele de co-expresie, accentul este pe **relațiile pereche gene–gene**, reprezentate ca muchii într-un graf, iar modulele sunt grupuri de gene puternic conectate între ele.

Clustering-ul răspunde la întrebarea „**ce gene/probe se grupează global împreună?**”,  
în timp ce rețeaua de co-expresie răspunde la „**cine este conectat cu cine și care gene par să funcționeze împreună în același modul?**”.  
Astfel, GCE-urile sunt mai potrivite pentru a identifica module funcționale locale, în timp ce clustering-ul oferă o vedere mai globală, dar mai rigidă, a structurii datelor.

