#rulare demo 1
Diagnosis,KMeans_Cluster
1,1
1,1
1,1
1,1
1,1
1,1
1,1
1,1
1,1
1,1
1,0
1,1
1,1
1,0
1,1
1,1
1,0
1,1
1,1
0,0
...

Am aplicat trei metode de clustering:
-hierarchical Clustering, oferă o imagine globală a similarității dintre probe, dar este mai greu de interpretat când sunt multe eșantioane 
-K-means (K=2) a separat destul de clar eșantioanele în două grupe corespunzătoare celor două diagnostice (malign și benign)
-DBSCAN a detectat câteva puncte anormale (outliers), dar nu a format clustere coerente, probabil din cauza distribuției compacte a datelor

Cea mai potrivită metodă a fost K-means, deoarece a reprodus cel mai bine separarea naturală dintre eșantioanele maligne și benigne și oferă o interpretare ușoară a centrelor de clustere.

Atât clustering-ul, cât și arborii filogenetici caută relații de similaritate între entități biologice, dar diferă în tipul de relații pe care le descriu. Pe scurt, clustering-ul descoperă asemănări statistice, iar arborii filogenetici reconstruiesc legături evolutive. Ambele vizualizează relații între entități, dar în contexte biologice diferite.
