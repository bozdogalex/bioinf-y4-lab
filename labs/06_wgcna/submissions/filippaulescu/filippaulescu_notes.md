## Notes — Gene Co-Expression Network

# Metrica folosita:  
Am utilizat corelatia Spearman, deoarece este potrivita pentru date RNA-Seq si surprinde relatii monotone chiar daca distributiile nu sunt liniare.

# Pragul ales:  
Pentru construirea retelei am aplicat un prag de |cor| >= 0.6, astfel incat in graf au ramas doar legaturile puternice dintre gene.

# Rezultate:  
Dupa preprocesare (log-transformare si filtrare dupa varianta) au ramas 36 de gene.  
Reteaua finala a avut 36 de noduri si 131 de muchii, iar algoritmul Louvain a identificat 6 module.  
Exemplu: gene ribozomale precum RPL13, RPL19, RPS23 si RPS29 au fost grupate in acelasi modul; la fel, gene imune precum NKG7 si CTSW au aparut impreuna.

# Cum difera o retea de co-expresie fata de clustering?  
Clustering-ul grupeaza genele doar in functie de asemanarea lor globala.  
O retea de co-expresie arata relatiile pereche dintre gene (muchii), iar modulele se formeaza pe baza conexiunilor din graf, nu doar pe similaritatea generala. Practic, clustering-ul creeaza grupuri "in masa", in timp ce retelele evidentiaza cine este conectat cu cine.
