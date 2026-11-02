# Lab 6 — Construirea rețelelor de co-expresie

> Soluțiile și PR-urile de evaluare rămân în repository-ul vostru **privat**, partajat cu `bozdogalex`. Trimiteți linkul PR și commit SHA în LMS. [Ghid de predare](../../docs/git-workflow.md).

După clustering păstrăm explicit relațiile pereche dintre gene, construim un graf și detectăm module. Lab 7 continuă cu **aceeași rețea salvată**, pentru vizualizare, hub-uri și interpretare. Practicul de aici folosește corelație, prag și Louvain; prezentarea compară această abordare cu WGCNA. Practicul nu implementează soft thresholding, TOM sau întregul pipeline WGCNA.

## Înainte de exercițiu

Aveți nevoie de operațiile pandas pentru indexare și transpunere, de corelația Pearson/Spearman și de reprezentarea unui graf neorientat prin adiacență. Fiecare linie din intrare este profilul unei gene pe aceleași probe. Un prag decide ce perechi devin muchii; ponderea unei muchii păstrate exprimă intensitatea asocierii. Louvain produce o partiție a genelor conectate, iar o genă fără muchii primește eticheta `-1`.

În tabelul de comparație, **densitatea** este proporția muchiilor existente din cele `n*(n-1)/2` perechi posibile, pentru toate genele păstrate, inclusiv izolate. **Cel mai mare modul** este numărul de gene din cea mai mare comunitate detectată; pentru o rețea fără muchii este zero. Numărul modulelor nu include eticheta `-1`.

## Datele folosite

Setul [data/sample/lab06](../../data/sample/lab06/) este disponibil offline: **120 de gene fictive × 48 de probe sintetice independente**, cu gene pe linii și probe pe coloane. `expression_matrix.csv` conține valori continue pe o scară de expresie deja pregătită, asemănătoare log2; nu sunt counts RNA-seq. Păstrăm valorile ca atare, fără încă o log-transformare sau normalizare. `sample_metadata.csv` descrie cele 24 de controale și 24 de probe tratate simulate.

Generatorul [generate_teaching_data.py](generate_teaching_data.py) folosește seed `20261101`. Are trei grupuri de câte 30 de gene și 30 de gene zgomot; un grup urmărește răspunsul opus altuia. `generation_groups.csv` descrie mecanismul simulării, nu modulele inferate. Valorile, parametrii și proveniența sunt documentate în `provenance.json`. Genele inventate nu au adnotări GO/KEGG și nu pot susține concluzii despre boli sau ținte terapeutice.

**Corectarea sursei vechi:** GSE115469 este un studiu single-cell al ficatului uman de la cinci donatori, nu un set de cancer mamar. Nu îl descărcăm pentru acest exercițiu. Mii de celule nu sunt mii de replici independente de donator; o analiză biologică a acelui set cere o alegere explicită a unității de analiză și preprocesare potrivită single-cell. [Înregistrarea GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469).

## Rulare în Codespaces

Din terminalul repository-ului vostru privat, înlocuiți `student` cu identificatorul vostru:

```bash
python labs/06_wgcna/demo01_corr_threshold.py
python labs/06_wgcna/ex01_gce_networks.py --handle student
```

Scriptul de bază este executabil. Îl rulăm pentru a verifica intrarea și a crea exportul folosit în Lab 7. Implementarea cerută studenților este apoi funcția de comparație din `ex02_compare_networks.py`; simpla rulare a baseline-ului nu încheie exercițiul. Opțiunile fiecărui script sunt disponibile prin `--help`.

Implicit: Spearman, relații pozitive `r >= 0.7`, muchii ponderate cu `r`, varianță minimă `0.05`, Louvain cu rezoluția `1.0` și seed `42`. Corelația între gene se calculează pe `df.T`, astfel încât rezultatul să aibă forma gene × gene. Diagonala adiacenței este zero.

- `--network positive`: păstrează relațiile pozitive care trec pragul.
- `--network unsigned`: prag pe `|r|`; răspunsurile opuse pot ajunge puternic conectate. Ponderea este `|r|`, iar semnul original rămâne în `correlation.csv` și `edges.csv`.
- `--binary`: atribuie ponderea 1 tuturor muchiilor păstrate.
- `--method pearson`: compară dependența liniară cu varianta Spearman bazată pe ranguri.
- `--input-scale prepared`: nu transformă din nou valorile. Opțiunea `normalized_counts` aplică `log2(x+1)` unor counts **deja normalizate**; nu face normalizarea între probe. Pentru date proprii, documentați pregătirea și unitatea biologică a coloanelor.

`positive` este un prag pe corelații pozitive, distinct de transformarea signed din WGCNA. Parametrul `--resolution` aparține obiectivului Louvain; valoarea sa schimbă scara comunităților.

O rulare identică poate fi repetată în același folder. Dacă schimbați datele, metadatele, parametrii sau versiunile bibliotecilor, baseline-ul refuză înlocuirea rezultatului existent. Folosiți un alt `--output-dir` pentru variante. `--overwrite` permite înlocuirea explicită; nu îl folosiți pentru comparații dacă vreți să păstrați baseline-ul pentru Lab 7.

## Comparațiile din laborator

1. Verificați calculul unei perechi din exemplul cu șase gene și identificați dacă ea trece pragul pozitiv `0.7`.
2. În `ex02_compare_networks.py`, implementați **numai `compare_networks`**. Funcția primește matricea de corelație deja calculată; nu încarcă datele și nu recalculează corelația în buclă. Pentru cele șase combinații dintre pragurile `0.5`, `0.7`, `0.9` și rețelele `positive`, `unsigned`, construiți adiacența, graful și partiția folosind funcțiile baseline-ului.
3. Întoarceți un DataFrame cu o linie pentru fiecare configurație și coloanele `Network`, `Threshold`, `Edges`, `Density`, `Isolates`, `Modules`, `LargestModule`. Folosiți ponderi și păstrați seed `42`, rezoluția `1.0`, aceeași metodă și aceleași gene. Excludeți `-1` din numărarea modulelor și a dimensiunilor lor. Driverul adaugă în tabel metoda, seed-ul, rezoluția și dimensiunea intrării.
4. Verificați rezultatul: fără self-loops, același număr de noduri și cel mult `n*(n-1)/2` muchii. La metodă fixă, numărul muchiilor poate doar să scadă sau să rămână egal când crește pragul. Pentru configurația implicită, numerele din tabel trebuie să coincidă cu `parameters.json` din baseline. Numărul modulelor nu are obligația să varieze monoton.
5. Explicați o schimbare concretă de conectivitate și diferența dintre rețeaua pozitivă și `|r|`. Comparați o genă dintr-un modul mare cu una rămasă izolată. Puteți consulta grupurile generatorului doar după inferență, pentru a explica efectele alegerii reprezentării.

Scaffold-ul conține un `NotImplementedError` explicit până când completați funcția; baseline-ul și demo-ul rulează independent. După implementare:

```bash
python labs/06_wgcna/ex02_compare_networks.py --handle student
```

Driverul calculează corelația o singură dată și salvează `comparison.csv` și `comparison_parameters.json` în `data/work/student/lab06/comparisons/`, fără să rescrie exportul de bază pentru Lab 7. Pentru toate cele șase variante se schimbă numai regula muchiilor. Dacă extindeți comparația la Pearson, calculați câte o matrice pentru fiecare metodă, în afara buclei pe praguri.

Pentru exportul complet al unei singure variante suplimentare puteți folosi un folder separat:

```bash
python labs/06_wgcna/ex01_gce_networks.py --handle student --threshold 0.5 --output-dir data/work/student/lab06/t05
python labs/06_wgcna/ex01_gce_networks.py --handle student --threshold 0.9 --output-dir data/work/student/lab06/t09
python labs/06_wgcna/ex01_gce_networks.py --handle student --network unsigned --output-dir data/work/student/lab06/unsigned
```

O genă izolată are `Module = -1`: este păstrată în export, dar nu formează un modul detectat. Modulele Louvain se calculează pe subgraful cu muchii. Etichetele numerice ale modulelor sunt identificatori, nu o ordonare biologică. Același seed reproduce aceeași inițializare în același mediu; păstrați versiunile și ordinea intrării pentru comparații reproductibile.

Pentru comparația cu mecanismul simulării puteți consulta grupurile generate; nu folosiți aceste etichete ca intrare în corelație sau Louvain. Folosiți varianta implicită pentru continuarea în Lab 7 și notați separat comparațiile.

## Exportul pentru Lab 7

Rularea implicită salvează în `data/work/student/lab06/`, folder local ignorat de Git:

| Fișier | Conținut |
| --- | --- |
| `processed_expression.csv` | Genele păstrate și valorile folosite în calcul |
| `correlation.csv` | Corelațiile cu semnul original |
| `adjacency.csv` | Matricea exactă a rețelei, cu aceleași gene pe linii și coloane |
| `edges.csv` | Perechi conectate, pondere și corelație originală |
| `modules.csv` | `Gene`, `Module`, `Isolate`, inclusiv genele izolate |
| `sample_metadata.csv` | Condiția probelor, în ordinea coloanelor expresiei |
| `parameters.json` | Prag, metodă, scară, seed, rezoluție, dimensiuni, versiuni și hash-ul intrării |

Lab 7 încarcă adiacența salvată; nu recalculează corelațiile cu parametri care ar putea produce altă rețea.

## Predare

În PR-ul privat includeți implementarea din `ex02_compare_networks.py`, celelalte modificări cerute la curs și `labs/06_wgcna/submissions/student/notes.md`. Explicația trebuie să conțină comenzile de reproducere, parametrii, hash-ul intrării, tabelul rezumativ al celor șase configurații, o pereche verificată și efectul schimbării pragului sau al folosirii `|r|`. Puteți include un eșantion de metadate de cel mult 10 linii, conform [politicii de date](../../docs/GDPR_and_DataPolicy.md).

Matricile complete, lista de muchii, mapping-ul modulelor și celelalte rezultate generate rămân în `data/work/student/lab06/`; nu le adăugați în Git, nici în PR-ul privat. Păstrați-le local pentru continuarea în Lab 7. Respectați și [politica privind asistența AI](../../docs/policies.md).

## Resurse

- [Fișa laboratorului](../../docs/lab_onepagers/06_WGCNA.md)
- [Continuarea — Lab 7](../07_network_viz/README.md)
- [Langfelder și Horvath, WGCNA (2008)](https://doi.org/10.1186/1471-2105-9-559)
- [FAQ WGCNA: date, eșantionare, preprocesare](https://edo98811.github.io/WGCNA_official_documentation/faq.html)
- [NetworkX: Louvain](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.community.louvain.louvain_communities.html)

Materialele de bază și corecturile pipeline-ului au fost pregătite cu asistență AI, sub îndrumarea cadrului didactic. Datele sunt o simulare documentată, nu rezultate experimentale.
