# Lab 6 — De la profiluri de expresie la module

Pornim de la un set **sintetic, pregătit**, cu 120 de gene fictive și 48 de probe independente. Genele sunt pe linii, probele pe coloane. Valorile sunt continue, pe o scară asemănătoare log2, nu counts RNA-seq. Nu aplicăm încă o log-transformare.

Cunoștințe necesare: indexare și transpunere pandas, corelații între profiluri, graf neorientat și matrice de adiacență. Densitatea este proporția perechilor de gene conectate; genele izolate rămân în numărul total de noduri. Un modul este o comunitate detectată, iar `-1` marchează o genă izolată.

## Traseul analizei

1. Citim matricea și eliminăm genele constante sau cu varianță sub prag.
2. Calculăm corelații între gene, pe toate cele 48 de probe.
3. Păstrăm muchiile care trec pragul: implicit `r >= 0.7`, cu ponderea `r`.
4. Construim graful neorientat, cu diagonala zero.
5. Detectăm module prin Louvain; păstrăm izolatele cu modul `-1`.
6. Exportăm rețeaua, expresia, modulele și parametrii pentru Lab 7.

Acesta este un pipeline **corelație → prag → Louvain**, distinct de WGCNA, unde apar transformarea prin putere, TOM și clusteringul ierarhic.

```bash
python labs/06_wgcna/demo01_corr_threshold.py
python labs/06_wgcna/ex01_gce_networks.py --handle student
```

Implementați `compare_networks` din `ex02_compare_networks.py`: folosiți aceeași corelație pentru cele șase combinații dintre pragurile 0.5/0.7/0.9 și rețelele pozitivă/unsigned. Calculați numărul muchiilor, densitatea, izolatele, modulele și dimensiunea celui mai mare modul. Nu includeți `-1` în module și nu recalculați corelația pentru fiecare prag. Driverul adaugă metoda, seed-ul și rezoluția și scrie numai comparația, fără să modifice baseline-ul.

```bash
python labs/06_wgcna/ex02_compare_networks.py --handle student
```

Scaffold-ul anunță explicit că funcția trebuie completată. Baseline-ul rămâne executabil. Pentru variante exportate complet folosiți `--output-dir` separat; schimbarea parametrilor în folderul existent cere `--overwrite` explicit.

Rularea implicită exportă `processed_expression.csv`, `correlation.csv`, `adjacency.csv`, `edges.csv`, `modules.csv`, `sample_metadata.csv` și `parameters.json` în `data/work/student/lab06/`, ignorat de Git. Păstrați aceste fișiere local pentru Lab 7; nu predați matricile sau listele complete.

În PR-ul privat includeți implementarea exercițiului și `labs/06_wgcna/submissions/student/notes.md`, cu comenzile, parametrii, hash-ul intrării și tabelul rezumativ al celor șase configurații. Opțional, un eșantion de metadate de cel mult 10 linii. [Politica de date](../GDPR_and_DataPolicy.md).

GSE115469, menționat în materialele vechi, este single-cell RNA-seq de ficat de la cinci donatori, nu cancer mamar; nu este intrarea acestui exercițiu. Identificatorii inventați din simulare nu se folosesc la GO/KEGG.

- [Instrucțiuni și parametri](../../labs/06_wgcna/README.md)
- [Date sintetice și proveniență](../../data/sample/lab06/)
- [Continuarea în Lab 7](../../labs/07_network_viz/README.md)
- [WGCNA, Langfelder și Horvath (2008)](https://doi.org/10.1186/1471-2105-9-559)

Material actualizat cu asistență AI, sub îndrumarea cadrului didactic.
