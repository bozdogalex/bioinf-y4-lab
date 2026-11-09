# Lab 7 — Interpretarea rețelelor de co-expresie

Continuați **varianta de bază salvată în Lab 6**: setul simulat de 120 de gene și 48 de probe, graful și modulele deja calculate. Folosiți [README-ul](README.md) pentru rulare și [prezentarea în română](../../docs/presentations/07_network_viz.pdf) pentru explicații. Înlocuiți `student` cu handle-ul folosit în Lab 6. Comenzile se rulează din rădăcina repository-ului privat.

## 1. Încărcarea rezultatului din Lab 6 — 2 puncte

```bash
python labs/07_network_viz/ex01_network_viz.py --handle student --seed 42 --output-dir data/work/student/lab07/seed42
```

Păstrați adiacența, modulele și parametrii aceleiași rulări. Pentru un director alternativ folosiți `--lab06-dir`. Dacă lipsește rezultatul de bază, îl regenerați prin `python labs/06_wgcna/ex01_gce_networks.py --handle student`.

În raport notați metrica, regula pentru semn, pragul, seed-ul Louvain, hash-ul intrării și numărul genelor, muchiilor, izolatelor și modulelor. Verificați concordanța cu `parameters.json`. Explicați ce reprezintă genele cu `Module = -1` și cum apar în vizualizare.

## 2. Compararea layout-urilor — 2 puncte

```bash
python labs/07_network_viz/ex01_network_viz.py --handle student --seed 7 --output-dir data/work/student/lab07/seed7
```

Comparați cele două fișiere `network.png` și cele două tabele `hubs.csv`. Precizați ce s-a schimbat și ce valori au rămas identice. Explicați de ce o genă care apare în centrul desenului poate avea un scor de conectivitate mai mic decât una aflată la margine. Păstrați graful fix pentru această comparație.

## 3. Analiza unui candidat hub — 4 puncte

Alegeți o genă conectată. Păstrați scripturile didactice nemodificate și scrieți `labs/07_network_viz/submissions/student/analyze_hub.py`. Programul propriu primește directorul exportului Lab 6, gena aleasă și un director local de ieșire. El trebuie să:

- încarce graful complet și modulele, inclusiv genele izolate;
- calculeze gradul genei, suma tuturor ponderilor și suma ponderilor către gene din același modul;
- extragă subgraful indus de genă și vecinii ei direcți și să îl salveze ca figură în `data/work/student/lab07/`;
- raporteze scorurile calculate pe graful complet, înaintea extragerii subgrafului.

Puteți reutiliza funcțiile didactice sau implementa citirea și calculele în programul propriu. Verificați scorurile față de `hubs.csv` dacă gena apare în lista exportată. În raport includeți scorurile, numărul vecinilor și interpretarea lor. Explicați diferența dintre un hub de co-expresie și o genă despre care avem dovezi experimentale că reglează alte gene.

Punctaj: citire și păstrarea genelor — 1p; cele trei scoruri — 1p; subgraf și export local — 1p; interpretare și reproducere — 1p.

## 4. Interpretarea modulului — 2 puncte

Formulați o concluzie susținută de simulare și precizați ce date suplimentare ar permite interpretarea biologică. Identificatorii `G001`–`G120` sunt fictivi. Un GO/KEGG autentic cere gene reale, o specie, identificatori recunoscuți și un univers eligibil documentat.

În exemplul ipotetic din prezentare, același modul are 30 de gene, dintre care 12 adnotate cu termenul T. Explicați efectul schimbării universului de la `(N=120, K=24)` la `(N=60, K=18)` asupra îmbogățirii. Precizați ce comparație schimbă universul și de ce testarea multor termeni cere o corecție.

În exemplul Diseasome, explicați cum o singură genă comună între trei boli produce trei muchii în proiecția pe boli și ce observații ar fi necesare pentru investigarea comorbidității.

## Extensie opțională — până la 1 punct bonus

```bash
python labs/07_network_viz/ex02_module_profiles.py --handle student
```

Scriptul complet executabil calculează eigengene-ul fiecărui modul și asocierea descriptivă cu tratamentul. Folosește `processed_expression.csv` și `sample_metadata.csv` din aceeași rulare Lab 6. Rezultatele sunt în `data/work/student/lab07/profiles/`.

Comparați modulul asociat pozitiv cu tratamentul, cel asociat negativ și modulul cu asociere slabă. Explicați cum poate un grup să aibă gene bine corelate între ele și o relație slabă cu tratamentul. Raportați proporția variației explicate de PC1 și convenția pentru semnul profilului. Bootstrap-ul, replicarea independentă și un GO/KEGG real sunt subiecte de extindere, în afara implementării obligatorii.

## Predare

Deschideți PR-ul **în repository-ul vostru privat**, partajat cu `bozdogalex`, și trimiteți în LMS linkul PR-ului și commit SHA, conform [ghidului de predare](../../docs/git-workflow.md).

Includeți în PR:

- `labs/07_network_viz/submissions/student/analyze_hub.py` și alte fișiere sursă proprii;
- `labs/07_network_viz/submissions/student/notes.md`, un raport concis cu comenzile, parametrii, hash-ul intrării și explicațiile cerute;
- atribuirea asistenței AI, dacă a fost folosită, conform [politicilor cursului](../../docs/policies.md).

Figurile, matricile, listele complete și tabelele generate rămân în `data/work/student/lab07/`, ignorat de Git. Predarea conține codul și explicațiile care permit reproducerea lor. Un eșantion opțional de metadate are cel mult 10 linii. Respectați [politica de date](../../docs/GDPR_and_DataPolicy.md).

Enunț actualizat cu asistență AI, sub îndrumarea cadrului didactic.


