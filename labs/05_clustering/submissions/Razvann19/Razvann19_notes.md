# Săptămâna 5 — Clustering în Bioinformatică

### Metoda cea mai potrivită
Pentru dataset-ul WDBC (cancer mamar), metoda **K-means** a oferit rezultate clare, separând probele maligne și benigne în două clustere bine delimitate după reducerea dimensionalității PCA.  
**DBSCAN** a detectat câteva puncte izolate, dar a considerat prea multe probe ca “noise”, iar **Hierarchical Clustering** a fost mai greu de interpretat vizual pentru acest tip de date cu multe dimensiuni.  
Prin urmare, **K-means (K=2)** s-a dovedit cea mai potrivită metodă pentru acest set de date.

### Comparație cu arborii filogenetici
Clustering-ul și arborii filogenetici sunt ambele metode de grupare, dar scopurile diferă:
- **Clustering-ul** caută **grupări pe baza similarității matematice** (ex. expresia genică, distanțe Euclidiene), fără o relație temporală sau evolutivă.  
- **Arborii filogenetici** reconstruiesc **relații de descendență evolutivă**, având o structură ierarhică direcționată.  
În concluzie, clustering-ul este util pentru **detectarea tiparelor funcționale sau moleculare** (ex. subtipuri de cancer), în timp ce filogenia explică **relațiile evolutive între gene sau organisme**.
