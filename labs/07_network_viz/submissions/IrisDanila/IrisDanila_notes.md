## Metodă de layout folosită

Am utilizat **spring layout** (algoritmul Fruchterman-Reingold) pentru vizualizarea rețelei de co-expresie genică. Această metodă poziționează nodurile astfel încât muchiile au lungimi aproximativ egale, iar nodurile conectate sunt atrase unele de altele, în timp ce nodurile neconectate se resping reciproc.

---

## Ce avantaje aduce vizualizarea față de analiza numerică din Lab 6? 

Analiza numerică din Lab 6 oferă precizie cantitativă (corelații, indici de modularitate), dar este dificil de interpretat intuitiv — o matrice de 1000×1000 sau o listă de module nu permit identificarea rapidă a structurii rețelei. 

**Vizualizarea aduce următoarele avantaje:**

1. **Identificare imediată a structurii modulare:** Modulele apar ca clustere vizuale de noduri colorate identic, făcând evident gradul de separare între grupuri funcționale.

2. **Recunoașterea genelor hub:** Genele cu grad mare (hub-uri) sunt vizibile ca noduri centrale în clustere, sugerând roluri regulatory importante — acestea pot fi ținte terapeutice în cancer sau alte boli.

3. **Validare vizuală a calității modulelor:** Dacă modulele sunt bine separate în spațiu, algoritmul de clustering a funcționat corect.  Module amestecate sugerează zgomot sau parametri neoptimali.

4. **Facilitarea comunicării interdisciplinare:** O figură este mult mai ușor de înțeles pentru biologi și clinicieni decât tabele numerice, accelerând interpretarea biologică și generarea de ipoteze.

5. **Contextualizare biomedicală (diseasome):** Vizualizarea permite suprapunerea informațiilor despre boli (ex: gene mutate în cancer) și identificarea modulelor asociate cu anumite patologii — util pentru drug repurposing și înțelegerea comorbidităților.

---

## Concluzie

Vizualizarea nu înlocuiește analiza numerică, ci o completează — oferă o perspectivă intuitivă care facili