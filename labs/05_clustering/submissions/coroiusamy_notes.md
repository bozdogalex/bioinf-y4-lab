# Laborator 5: Clustering - Analiză și Reflecții

## 1. Metoda cea mai potrivită pentru datele analizate

Pentru setul de date WDBC (cancer mamar), cea mai potrivită metodă de clustering s-a dovedit a fi **K-means (cu K=2)**.

**Justificare:**

* **Potrivirea cu problema:** Știam din contextul problemei (diagnostic Malign vs. Benign) că ne așteptăm să găsim **două** grupuri principale. K-means este optimizat exact pentru acest scenariu, în care numărul de clustere (K) este cunoscut apriori.
* **Rezultate clare:** Vizualizarea PCA a arătat o separare foarte clară și distinctă a celor două clustere găsite de K-means, sugerând că algoritmul a capturat cu succes structura de bază a datelor.
* **Confirmare:** Clustering-ul Ierarhic (Hierarchical Clustering) a susținut această concluzie, dendrograma arătând clar două ramuri principale care se unesc la o distanță mare, confirmând existența a două super-grupuri naturale.
* **Limitarea DBSCAN:** În contrast, DBSCAN a eșuat să găsească o structură relevantă cu parametrii dați (`eps=1.5, min_samples=5`), clasificând majoritatea punctelor drept zgomot (eticheta -1). Acest lucru evidențiază dificultatea alegerii parametrilor corecți pentru DBSCAN pe un set de date necunoscut, făcându-l mai puțin potrivit pentru această analiză exploratorie inițială.

---

## 2. Reflecție: Clustering vs. Arbori Filogenetici

Deși ambele metode grupează datele pe baza similarității și pot produce vizualizări arborescente (o dendrogramă seamănă cu un arbore filogenetic), scopul și interpretarea lor sunt fundamental diferite în biologie.

* **Arborii Filogenetici (Săpt. 4):**
    * **Scop:** Descoperirea **relațiilor evolutive** și a **istoriei ancestrale**.
    * **Întrebare:** "Cine este strămoșul comun?" sau "Cât de recent au divergat două specii/gene?".
    * **Axe:** Axele reprezintă **timpul evolutiv** sau **distanța genetică** (nr. de mutații).
    * **Aplicație:** Se folosesc pe secvențe ADN/proteine pentru a înțelege cum au evoluat organismele.

* **Clustering (Săpt. 5):**
    * **Scop:** Descoperirea **tiparelor (patterns)** și a **grupurilor de similaritate** în datele actuale.
    * **Întrebare:** "Ce pacienți/gene se comportă similar *acum*?" sau "Există sub-tipuri de cancer în aceste date?".
    * **Axe:** Axele unei dendrograme reprezintă o **distanță matematică** (ex: euclidiană), nu neapărat timp.
    * **Aplicație:** Se folosesc pe date de expresie genică, profile metabolice etc., pentru a găsi grupuri funcționale (ex: gene co-exprimate) sau sub-tipuri de boli (ex: clustere de pacienți).

**Pe scurt:** Filogenetica ne spune despre **înrudire și istorie (ancestry)**, în timp ce clustering-ul ne spune despre **similaritate și funcție (patterns)**.