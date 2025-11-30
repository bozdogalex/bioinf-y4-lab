**Metoda folosită**

Pentru vizualizarea rețelei am folosit layout-ul spring_layout, care aranjează automat nodurile astfel încât genele conectate să stea aproape unele de altele.Nodurile au fost colorate după modulul din care fac parte și mărite dacă aveau multe conexiuni, așa încât genele cele mai importante să iasă ușor în evidență în rețea.

**Ce avantaje aduce vizualizarea față de analiza numerică din Lab 6?**  

Vizualizarea este mult mai intuitivă decât analiza numerică din laboratorul 6, deoarece putem vedea imediat cum se organizează genele, cât de dense sunt modulele și ce gene ies în evidență ca hub-uri.
De exemplu, din vizualizare se vede imediat că unele module stau strânse într-un singur grup, genele imune (CST7, CTSW, NKG7) sunt apropiate și bine conectate între ele, iar gene precum APCS sau APOH ies în evidență prin numărul mare de legături din modulul lor - aspecte care nu sunt atât d evidente doar din tabele sau valori numerice.

**Bonus**

- Vizualizare în Cytoscape
După ce am construit rețeaua și am generat fișierul cu muchii (edges_eeeevn0.csv), l-am încărcat în Cytoscape, unde nodurile au fost aranjate mult mai clar, iar modulele au fost separate vizibil, ceea ce a făcut rețeaua mult mai ușor de înțeles.

- Comparație cu NetworkX
Față de imaginea statică generată de NetworkX, în Cytoscape am putut mări și explora zonele rețelei, iar grupuri precum genele ribozomale și cele imune se văd mult mai clar și sunt mult mai ușor de recunoscut.