**Demo01:**

Am rulat demo01, am pus rezultatul aici ca fișierul `clusters_demo.csv`

**Task01:**

Am completat `ex01_clustering.py` și am obțiunut cele 3 documente .png și cel de tip .csv.

Dintre metodele testate (Hierarchical, K-means și DBSCAN), *K-means* s-a dovedit cea mai potrivită pentru acest set de date.  
Datele din WDBC sunt destul de bine separate în două grupe (maligne și benigne), iar K-means, având un număr fix de clustere (K=2), a reușit să le distingă clar în planul PCA.  
Hierarchical clustering oferă o privire interesantă asupra relațiilor dintre mostre, dar tinde să devină greu de interpretat pentru un număr mare de eșantioane.  
DBSCAN a detectat unele puncte ca “zgomot” (noise), însă nu a reușit să grupeze complet datele, semn că nu este metoda ideală pentru acest tip de distribuție compactă.

*Concluzie:* K-means este cea mai stabilă și interpretabilă metodă pentru acest dataset numeric standardizat, oferind separații coerente între clasele biologice.

Atât *clustering-ul*, cât și *arborii filogenetici* caută să grupeze entități în funcție de asemănări, însă scopurile lor diferă:
- Clustering-ul (K-means, DBSCAN etc.) organizează datele doar pe baza *similarității* între trăsături numerice (ex: dimensiunea celulelor, textura, forma).  
- Arborii filogenetici, în schimb, reprezintă *relații evolutive* și implică o direcție istorică, sugerând o „descendență comună”.  
- Dendrograma din hierarchical clustering seamănă vizual cu un arbore filogenetic, dar *nu exprimă evoluția*, ci doar gruparea pe baza distanțelor matematice.

Astfel, clustering-ul este o metodă de *organizare a datelor*, pe când arborii filogenetici sunt o cale de *a interpreta biologic relațiile dintre organisme*.
