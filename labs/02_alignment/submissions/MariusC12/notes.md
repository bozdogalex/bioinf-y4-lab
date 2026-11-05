# Notes – Lab 2 Sequence Alignment

## Task 1 – Distanțe p

Am calculat distanțele p pentru toate cele trei perechi de secvențe folosind alinieri globale.

Rezultatele obținute:

| Pereche | Mismatch-uri | Poziții comparate | Distanță p |
|---|---:|---:|---:|
| Human – Mouse | 287 | 1749 | 0.1641 |
| Human – Zebrafish | 836 | 2195 | 0.3809 |
| Mouse – Zebrafish | 648 | 1755 | 0.3692 |

Cea mai mică distanță a fost obținută între Homo sapiens și Mus musculus, ceea ce indică o asemănare mai mare între aceste secvențe. Cea mai mare distanță a fost între Homo sapiens și Danio rerio.

La calcul au fost excluse pozițiile cu goluri și nucleotide ambigue. Numărul pozițiilor comparate diferă între perechi.

Distanța p reprezintă proporția diferențelor observate. Aceasta nu este o distanță evolutivă corectată și nu demonstrează singură o relație filogenetică.

## Task 2 – Aliniere globală și locală

Pentru comparație am ales două fragmente CDS de câte 100 nucleotide din secvențele TP53 umană și de șoarece.

Accession-uri:
- Human: NM_000546.6
- Mouse: NM_011640.3

Coordonatele CDS:
- Human: [142:1324]
- Mouse: [157:1330]

Blocul aliniat din care au fost selectate fragmentele:
- Human CDS: [270:872]
- Mouse CDS: [261:863]

Coordonatele fragmentelor în CDS:
- Human: [521:621]
- Mouse: [512:612]

Coordonatele corespunzătoare în secvențele FASTA:
- Human: [663:763]
- Mouse: [669:769]

Am utilizat scorurile match = +1, mismatch = -1 și gap = -2.

### Rezultate

| Algoritm | Scor propriu | Scor Biopython |
|---|---:|---:|
| Needleman–Wunsch | 70 | 70 |
| Smith–Waterman | 72 | 72 |

Rezultatele implementărilor proprii coincid cu cele obținute folosind Biopython.

Needleman–Wunsch realizează o aliniere globală, luând în considerare întreaga lungime a celor două fragmente. Smith–Waterman realizează o aliniere locală și identifică regiunea cu cel mai bun scor.

În cazul nostru, alinierea locală are lungimea de 98 de poziții, în timp ce fragmentele inițiale au câte 100 de nucleotide.

Alinierea globală este utilă pentru compararea secvențelor pe toată lungimea, iar cea locală pentru identificarea unor regiuni similare în secvențe care pot conține și porțiuni diferite.


# Task 3 – aliniere multiplă (MSA)


## Task 3 – Aliniere multiplă (MSA)

Am realizat alinierea multiplă a celor trei secvențe de ARN mesager folosind Clustal Omega versiunea 1.2.4. Secvențele provin de la Homo sapiens, Mus musculus și Danio rerio.

Alinierea obținută conține 2628 de coloane. Am identificat un fragment conservat între coloanele 864 și 924. În acest fragment există 40 de poziții în care toate cele trei secvențe au aceeași nucleotidă.

Conservarea acestor poziții arată că există regiuni asemănătoare ale genei TP53 între cele trei organisme. Acest lucru poate sugera importanța biologică a regiunii respective, însă alinierea singură nu demonstrează funcția exactă a acesteia.

Secvența de la Danio rerio adaugă o perspectivă asupra conservării TP53 la o specie mai îndepărtată evolutiv față de om și șoarece. Compararea celor trei secvențe permite observarea asemănărilor și diferențelor dintre organisme.

Fragmentul MSA selectat a fost salvat în fișierul `task3_msa_excerpt.fasta`, iar alinierea completă în `task3_msa.fasta`.