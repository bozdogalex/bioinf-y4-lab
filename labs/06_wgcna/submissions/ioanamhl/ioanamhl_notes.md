-am rulat fisierul ex07_gce_networks.py cu fisierul expression_matrix2 ce contine primele 2000 de instante din matrice de expresie corespunzatoare genei GSE115469

Setari utilizate:
-metrica de corelație: Spearman  
-prag de corelație: 0.5  
-filtrare gene: prag = 0.6  
-tip rețea: neorientata  
-algoritm de detecție a modulelor: greedy modularity communities (NetworkX)  

Reflecție: rețele de co-expresie vs. clustering clasic

Clustering-ul clasic grupează genele pe baza asemănării profilelor de expresie, fără a descrie relațiile directe dintre ele. Rețelele de co-expresie modelează aceste relații folosind corelații între gene, permițând identificarea modulelor de gene co-reglate. Rezultatele depind însă de pragurile alese și pot fi sensibile la zgomotul din date.
