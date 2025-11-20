# Saptamana 5: Clustering in Bioinformatica - Analiza WDBC (Cancer Mamar)

## Rezumat

Au fost aplicate trei metode de clustering (Hierarchical, K-Means, DBSCAN) pe datele standardizate WDBC, care contin 30 de caracteristici masurate pe nucleii celulari. Scopul a fost separarea automata a esantioanelor in grupuri, fara a folosi eticheta reala de diagnostic.

### 1. K-Means (K=2)
Rezultatul a aratat o separare vizuala clara intre doua grupuri. Prima componenta principala (PC1) diferentiaza cel mai bine cele doua tipuri de celule, probabil Benigne si Maligne.

### 2. Hierarchical (Average Linkage)
Dendrograma arata doua grupuri principale daca se aplica un prag la distanta de aproximativ 12.5. Acestea corespund in mare cu separarea observata in K-Means.

### 3. DBSCAN (eps=1.5, min_samples=5)
DBSCAN a detectat un singur cluster mare si cateva puncte izolate. Valoarea eps a fost prea mare, ceea ce a facut ca cele doua grupuri sa fie unite. Parametrii trebuie ajustati mai fin pentru a evidentia doua clustere.

---

## Metoda cea mai potrivita

Cea mai potrivita metoda pentru aceste date este **K-Means (K=2)**, deoarece:
1. Ofera cea mai buna separare vizuala.
2. Se potriveste cu informatia biologica cunoscuta (Benign vs. Malign).
3. Este mai simpla si mai interpretabila decat celelalte metode.

---

## Clustering vs. Arbori Filogenetici

Clustering-ul grupeaza esantioanele in functie de similaritatea masuratorilor (ex. expresie genica, proprietati celulare). Arborii filogenetici cauta relatii evolutive intre secvente (ADN, ARN, proteine).

Clustering-ul este folosit pentru identificarea de subtipuri de boala sau gruparea pacientilor cu profiluri similare. Arborii filogenetici sunt folositi pentru a urmari originea comuna a speciilor sau genelor.

---

### Concluzie

Clustering-ul descrie asemanari functionale prezente intre esantioane, iar arborii filogenetici descriu originea lor evolutiva. Cele doua metode se completeaza in analiza bioinformatica.
