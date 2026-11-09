# Lab 7 — Interpretarea rețelei construite în Lab 6

> Soluțiile și PR-urile de evaluare rămân în repository-ul vostru **privat**, partajat cu `bozdogalex`. Trimiteți linkul PR și commit SHA în LMS. [Ghid de predare](../../docs/git-workflow.md).

Continuăm analiza din [Lab 6](../06_wgcna/README.md) folosind **același graf și aceleași module**. Vizualizarea schimbă pozițiile nodurilor în pagină, fără să recalculeze corelațiile sau muchiile.

## Rulare în Codespaces

După rularea variantei de bază din Lab 6, din rădăcina repository-ului privat:

```bash
python labs/07_network_viz/ex01_network_viz.py --handle student
```

Implicit, scriptul citește `adjacency.csv`, `modules.csv` și `parameters.json` din `data/work/student/lab06/`, unde au rămas salvate local. Dacă ați salvat într-un alt folder, specificați-l:

```bash
python labs/07_network_viz/ex01_network_viz.py --handle student --lab06-dir data/work/student/lab06/t05
```

Loaderul verifică simetria, diagonala zero, concordanța genelor și numerele de noduri/muchii/izolate/module față de `parameters.json`. Astfel poate detecta fișiere amestecate din rulări diferite. Genele izolate rămân în graf, au `Module = -1` și sunt desenate gri într-un inel exterior. Acest inel este o convenție de afișare, fără semnificație biologică.

## Ce analizăm

- Comparăm aspectul rețelei când schimbăm seed-ul layout-ului, păstrând graful fix. Distanța din desen nu este o măsurare a distanței dintre profilurile de expresie.
- Nodurile sunt colorate după modul și dimensionate după grad. Hub-urile exportate sunt ordonate după **conectivitatea ponderată în interiorul modulului**; gradul și suma ponderilor globale sunt raportate separat.
- Investigăm o genă bine conectată și vecinii ei. Conectivitatea măsoară poziția în rețeaua de co-expresie aleasă, nu dovedește reglare directă sau cauzalitate.
- Discutăm eigengene-ul unui modul ca prima componentă principală a profilurilor genelor sale și folosirea lui în comparații cu metadatele probelor.
- Discutăm validarea prin reeșantionarea probelor și prin date independente. O împărțire stabilă și o interpretare biologică plauzibilă răspund la întrebări diferite.

Acesta este un baseline executabil pentru explorare. Parametrii sunt disponibili prin `--help`. Păstrați scripturile didactice nemodificate și realizați analiza proprie în `submissions/<handle>/`, conform [enunțului și punctajului](assignment.md).

## Practicul din prezentarea nouă

Partea obligatorie folosește același graf pentru două layout-uri, apoi analizează o genă hub și vecinătatea ei. Rulați în directoare separate pentru a păstra ambele rezultate:

```bash
python labs/07_network_viz/ex01_network_viz.py --handle student --seed 42 --output-dir data/work/student/lab07/seed42
python labs/07_network_viz/ex01_network_viz.py --handle student --seed 7 --output-dir data/work/student/lab07/seed7
```

Comparați `network.png` și `hubs.csv`. Seed-ul de aici controlează numai pozițiile. Mărimea nodului reprezintă gradul, iar clasamentul hub-urilor folosește conectivitatea ponderată intramodulară. Calculați scorurile genei pe graful complet, înainte de extragerea subgrafului pentru afișare. Programul propriu și interpretarea se cer în [assignment.md](assignment.md).

## Extensie executabilă: eigengene și tratament

```bash
python labs/07_network_viz/ex02_module_profiles.py --handle student
```

Această extensie opțională citește și `processed_expression.csv` și `sample_metadata.csv` din rularea Lab 6. Fiecare genă este standardizată pe probe cu deviație standard de eșantion; SVD produce PC1 pentru fiecare modul conectat. Scorurile sunt orientate spre profilul mediu al genelor și standardizate la deviație standard unu. Genele izolate sunt păstrate în rețeaua de bază și excluse din calculul profilurilor de modul.

Metadatele sunt aliniate după identificatorii probelor. Exemplul așteaptă etichetele `control` și `tratament` din datele cursului. Corelația Pearson dintre profil și condiția codificată 0/1 este descriptivă, calculată pe probele folosite la descoperirea modulelor.

În `data/work/student/lab07/profiles/` se salvează `module_profiles.csv`, `module_summary.csv`, `module_profiles.png` și `profile_parameters.json`. Rezumatul include numărul genelor, proporția variației explicate de PC1 și corelația cu tratamentul. Opțiunile `--lab06-dir` și `--output-dir` permit selectarea altor directoare. Fișierele de expresie și metadate trebuie să provină din aceeași analiză.

Analiza funcțională, bootstrap-ul și replicarea în date independente sunt explicate în prezentare. Implementarea lor completă rămâne o extindere; exemplul obligatoriu de îmbogățire este un calcul pe adnotări ipotetice.

## Interpretarea funcțională și Diseasome

Setul implicit folosește **identificatori de gene inventați**. Nu îi trimiteți la GO/KEGG, g:Profiler sau baze de boli. Pentru simulare putem compara modulele inferate cu mecanismul generatorului, nu putem face enrichment biologic. O demonstrație de enrichment folosește separat o listă de gene reale, o specie și un background eligibil explicit; când testăm mai mulți termeni raportăm corecția pentru testări multiple. Rezultatele unei astfel de demonstrații nu se atribuie simulării.

Diseasome este o schimbare a reprezentării: de exemplu, bolile devin noduri și sunt conectate prin gene asociate comune. Muchiile sale au altă definiție decât corelațiile dintre expresiile genelor. Rețeaua bolilor nu rezultă automat din modulele de co-expresie; asocierea la boli cere surse externe și o regulă explicită de conectare. Nici o legătură prin gene comune nu stabilește singură o comorbiditate observată clinic.

[Prezentarea în română](../../docs/presentations/07_network_viz.pdf) continuă Lab 6 cu vizualizarea, conectivitatea, profilurile modulelor, interpretarea funcțională și modelul Diseasome. [Fișa scurtă](../../docs/lab_onepagers/07_network_viz.md) și [enunțul](assignment.md) folosesc aceleași intrări și convenții.

## Export și predare

Rezultatele se salvează în `data/work/student/lab07/`, folder local ignorat de Git:

- `network.png` — graful cu modulele, hub-urile etichetate și izolatele vizibile;
- `hubs.csv` — `Gene`, `Module`, `Degree`, `Strength`, `IntramodularStrength`;
- `visualization_parameters.json` — seed și convenții de desen, împreună cu parametrii Lab 6;

În PR-ul privat includeți programul propriu `labs/07_network_viz/submissions/student/analyze_hub.py` și `labs/07_network_viz/submissions/student/notes.md`, cu comenzile de reproducere, parametrii rețelei și layout-ului, hash-ul intrării, un rezumat al dimensiunilor rețelei și explicațiile din enunț. Trimiteți linkul PR și commit SHA în LMS. Opțional, includeți un eșantion de metadate de cel mult 10 linii. Matricile, listele complete, PNG-ul și celelalte rezultate generate rămân local; nu le adăugați în PR. Urmați [politica de date](../../docs/GDPR_and_DataPolicy.md).

Un export suplimentar în Cytoscape sau Gephi poate folosi `edges.csv` și `modules.csv` din Lab 6; în import se păstrează separat și genele izolate. Aceste fișiere rămân tot în zona locală de lucru.

## Resurse

- [NetworkX — desenarea grafurilor](https://networkx.org/documentation/stable/reference/drawing.html)
- [Cytoscape](https://cytoscape.org/)
- [g:Profiler](https://biit.cs.ut.ee/gprofiler/)
- [Goh et al., The human disease network, PNAS (2007)](https://doi.org/10.1073/pnas.0701361104)
- [Următorul laborator — Machine Learning](../08_ML_flower/README.md)

Corecturile handoff-ului și materialele de bază au fost pregătite cu asistență AI, sub îndrumarea cadrului didactic.
