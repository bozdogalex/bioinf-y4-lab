# Saptamana 5: Clustering in Bioinformatica - Analiza WDBC (Cancer Mamar)

## Rezumatul Analizei

Am aplicat trei metode de clustering (Hierarchical, K-means, DBSCAN) pe datele standardizate de cancer mamar (WDBC), care contin 30 de caracteristici masurate pe nucleii celulari. Obiectivul a fost separarea automata a esantioanelor in grupuri distincte, fara a folosi eticheta de diagnostic reala.

### 1. K-Means Clustering (K=2)
**Rezultat (Vizualizare PCA):** Metoda K-Means, fortata la $K=2$, a produs o **separare vizuala foarte buna** a esantioanelor. Cele doua clustere sunt bine definite in spatiul redus de PCA:
* Clusterul din stanga (galben/verde) se afla predominant in zona PC1 negativa.
* Clusterul din dreapta (violet/negru) se afla predominant in zona PC1 pozitiva.

Acest rezultat sugereaza ca prima componenta principala (PC1) capteaza cea mai mare parte a variatiei care distinge, cel mai probabil, celulele Benigne de cele Maligne.

### 2. Hierarchical Clustering (Average Linkage)
**Rezultat (Dendrograma):** Dendrograma arata o structura ierarhica clara. Daca s-ar trasa o linie orizontala (cut-off) la o distanta de aproximativ **12.5**, s-ar obtine **doua clustere principale** distincte. Aceste doua clustere corespund, cel mai probabil, separarii Benign/Malign observata in K-Means. O distanta mai mica (ex: 5.0) ar genera un numar mare de clustere mai mici si mai putin informative.

### 3. DBSCAN Clustering ($\epsilon=1.5, min\_samples=5$)
**Rezultat (Vizualizare PCA):** Cu hiperparametrii $\epsilon=1.5$ si $min\_samples=5$ pe datele standardizate si vizualizate prin PCA, DBSCAN a identificat **un singur cluster mare** si cateva puncte de zgomot izolate, desi acestea nu sunt colorate separat pe plot.
* **Interpretare:** Aceasta indica faptul ca, desi exista o separare vizibila pe PCA (ca in K-Means), densitatea datelor este continua. Distanta $\epsilon$ folosita a fost prea mare, permitand algoritmului sa "uneasca" cele doua grupuri, care sunt de fapt separate prin densitate scazuta, nu prin goluri mari. DBSCAN necesita o reglare fina a parametrilor pentru a descoperi structura de doua clustere.

---

## 🥇 Metoda cea mai potrivita pentru WDBC

Metoda pe care am considerat-o cea mai potrivita pentru a descoperi structura binara din datele WDBC este **K-means Clustering (K=2)**.

**Motivatie:**
1.  **Acuratete vizuala:** K-Means a realizat cea mai curata si logica separare vizuala in cele doua grupuri asteptate.
2.  **Cunoastere a priori:** Stim ca datele provin din doua clase biologice distincte (Benign si Malign). K-Means functioneaza excelent atunci cand numarul de clustere ($K=2$) este cunoscut sau poate fi presupus in mod rezonabil, reusind sa optimizeze centrele clusterelor pentru a maximiza separarea.
3.  **Dendrograma** a confirmat existent
