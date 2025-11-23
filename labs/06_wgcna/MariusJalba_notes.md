Cum diferă o rețea de co-expresie față de clustering-ul clasic (Lab 5)?

Clusterring-ul clasic grupeaza genele pe baza asemanarii globale dintre profilele lor de expresie. Metodele precum K-means sau clustering ierarhic analizeaza fiecare gena ca pe un vector si incearca sa identifice grupuri cu tipar similar la nivel general, acest tip de analiza este util pentru o segmentare rapida, dar nu capteaza complexitatea relatiilor dintre gene.

Reteaua de co-expresie porneste de la relatiile pereche dintre gene, construind un graf in care nodurile reprezinta gene, iar muchiile corelatii puternice de expresie. Aplicand asupra grafului algoritmi de detectare a comunitatilor, care identifica module dense de interactiuni. Aceasta metoda evidentiaza grupuri de gene care functioneaza impreuna si gene care au rol potential central in procese biologice.


Am folosit corelatia Spearman, deoarece surprinde relatiile monotone dintre gene si este robusta la distrubutii cu forma neregulata si valori extreme din RNA-Seq.
Matricea de adiacenta a fost construita folosind pragul 0.5 pastrand doar conexiunile cu co-expresie puternica si eliminand legaturile slabe
Filtrarea genelor cu varianta scazuta si aplicarea pragului au redus zgomotul si au generat o retea mai precisa. Aplicarea Louvain a identificat mai multe noduri distincte, reflectand grupuri de gene cu relatii stranse de co-expresie