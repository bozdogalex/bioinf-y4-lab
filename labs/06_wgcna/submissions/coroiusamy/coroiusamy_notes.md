
## Parametrii folosiți

Pentru construcția rețelei de co-expresie din `ex01_gce_networks.py`, am folosit următorii parametri principali:

* **Metrica de Corelație:** `spearman`
* **Prag de Adiacență:** `0.6` (aplicat pe valoarea absolută a corelației)
* **Prag de Varianță:** `0.1` (pentru a filtra genele înainte de eșantionare)
* **Eșantionare:** Am folosit un eșantion aleatoriu de **500 de gene** (dintr-un subset inițial de 10.000) pentru a permite rularea scriptului în limitele de memorie ale mediului.

---

## Reflecție: Rețele de Co-expresie vs. Clustering Clasic

Întrebarea: *Cum diferă o rețea de co-expresie față de clustering-ul clasic (Lab 5)?*

Deși ambele metode au ca scop gruparea genelor care se comportă similar, ele o fac în moduri fundamental diferite și răspund la întrebări diferite:

### 1. Nivelul de analiză: Global vs. Pereche

* **Clustering-ul Clasic (ex: K-Means, Clustering Ierarhic)** grupează genele pe baza similarității **globale** a profilurilor lor de expresie. O genă este comparată cu un "centru" al unui cluster (K-Means) sau cu toate celelalte gene (HCA) folosind o metrică de distanță (ex: Euclideană). Rezultatul este o atribuire "forțată": fiecare genă aparține unui singur cluster.
    * **Analogie:** Sortarea elevilor în "case" (ex: "sportivi", "artiști", "matematicieni") pe baza performanței lor generale.

* **Rețelele de Co-expresie (GCEs)** se concentrează pe **relațiile pereche (gene-la-gene)**. Ele calculează o corelație pentru *fiecare pereche de gene* și construiesc o rețea doar din perechile care depășesc un prag.
    * **Analogie:** Crearea unei rețele sociale a elevilor, unde o conexiune există doar dacă doi elevi vorbesc *direct* unul cu celălalt (sunt "corelați").

### 2. Rezultatul final

* **Clustering-ul** produce "containere" sau "cutii" (clusterele). O genă este fie înăuntru, fie în afară. Nu oferă informații despre *de ce* sau *cum* genele *din interiorul* aceluiași cluster sunt legate între ele.

* **Rețelele** produc o "hartă" a conexiunilor. Această hartă ne permite să:
    * **Detectăm Module:** Acestea sunt comunități *dens conectate* în cadrul rețelei (grupuri de prieteni apropiați în rețeaua socială).
    * **Identificăm Topologia:** Putem vedea structura rețelei. Putem găsi **gene "hub"** (gene centrale, extrem de conectate, ca "liderii" grupurilor de prieteni) sau gene "pod" care conectează module diferite.

**Pe scurt:** Clustering-ul ne spune **"ce grupuri există"** pe baza unei vederi de ansamblu, în timp ce rețelele de co-expresie ne spun **"cine este conectat cu cine"** și ne permit să deducem grupurile (modulele) pe baza densității acestor conexiuni.