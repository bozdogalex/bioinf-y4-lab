**Cea mai potrivita metoda**

K-Means (cu K=2) este cea mai potrivită.
Știm de la început că datele au exact două categorii (Malign și Benign). K-Means este conceput exact pentru a împărți datele în K grupuri clare, unde K este cunoscut.
    -> Clustering-ul Ierarhic este util pentru explorare (ca în dendrogramă), dar K-Means oferă o soluție finală mai directă.
    -> DBSCAN este nepotrivit aici. Este greu de reglat și este conceput pentru a găsi clustere de forme ciudate și a izola zgomotul (outlierii), ceea ce nu este cazul aici.


**Clustering vs. Arbori Filogenetici**

Deși dendrograma și un arbore filogenetic arată asemănător, ele măsoară lucruri diferite:
    -> Clustering-ul (din exercițiu): Grupează după similaritatea actuală a trăsăturilor.
        Exemplu: Cât de asemănătoare sunt două tumori pe baza măsurătorilor (rază, textură etc.)?
    -> Arborii Filogenetici: Grupează după istorie evolutivă și strămoși comuni.
        Exemplu: Cât de înrudite genetic sunt două specii (cât de recent au avut un strămoș comun)?
Pe scurt, clustering-ul arată cât de asemănătoare sunt lucrurile acum, în timp ce filogenetica arată relația lor de rudenie de-a lungul timpului.