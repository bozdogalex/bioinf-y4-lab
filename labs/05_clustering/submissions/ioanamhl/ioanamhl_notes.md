1. Dintre metodele testate K-means (K=2) a oferit cele mai clare rezultate pentru dataset-ul WDBC. 
Motivele:
- datele sunt standardizate și formează două grupuri relativ compacte (benign/malign);
- K=2 are sens biologic deoarece diagnosticul real are două clase;
- vizualizarea cu PCA arată clustere bine separate;
- metoda este stabilă și ușor de interpretat.

2. Cum se compară clustering-ul cu arborii filogenetici în descoperirea relațiilor biologice?
 Clustering:
- împarte datele în grupuri pe baza unei măsuri de similaritate;
- descoperă structuri ascunse fără presupuneri despre evoluție;
- este potrivit pentru expresie genică, date clinice, fenotipuri, transcriptomică.
Filogenetica:
- construiește un arbore care modelează relații ancestrale;
- are o structură direcțională și temporală (strămoș → descendent);
- aplicabilă la secvențe ADN/ARN, variante genomice, evoluție virală/bacteriană.

