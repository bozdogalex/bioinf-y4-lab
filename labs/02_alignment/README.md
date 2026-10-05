# Laboratorul 2 Alinierea secvențelor

Prerechizite: noțiunile introductive de genetică și genomică discutate la curs, formatele FASTA/GenBank și operațiile din Lab 1. Învățăm global, local, scoruri, gap-uri și programare dinamică, apoi interpretăm comparații biologice.

## Demonstrații

Din rădăcina repository-ului, în terminalul Codespaces:

```bash
python labs/02_alignment/demo01_pairwise.py --fasta data/sample/toy_alignment.fasta --k 10
python labs/02_alignment/demo01_pairwise.py --fasta data/sample/tp53_dna_multi.fasta --k 60
python labs/02_alignment/demo02_distance_matrix.py --fasta data/sample/tp53_dna_multi.fasta
```

Perechea artificială produce global −2 și local 4 cu scorurile demonstrației (+1, −1, gap −1). Demo-ul pairwise ia doar primele două înregistrări și primele `k` baze. Aceste prefixe nu sunt automat regiuni codante sau domenii.

Demo-ul de distanțe aliniază global fiecare pereche înainte de calcul. Exclude coloanele cu gap-uri și baze ambigue, raportează numărul pozițiilor comparate și `NA` când acesta este zero. Prin urmare, o pereche care diferă numai printr-un indel poate avea p-distance zero pe bazele comparate. Citiți alinierea și numitorul împreună cu rezultatul. Nu este o distanță evolutivă corectată și nu construim o filogenie doar din această demonstrație.

Pentru un FASTA deja aliniat, folosiți `--aligned`. Rândurile trebuie să aibă aceeași lungime; această condiție singură nu dovedește că secvențele au fost aliniate. Politica este pairwise deletion, deci numitorul poate diferi între perechi.

## Exerciții

Copiați `ex01_global_nw.py` și `ex02_local_sw.py` în `submissions/<handle>/`. Completați inițializarea și recurența. Nu modificați scheletele originale. Începeți cu `data/sample/toy_alignment.fasta`, apoi cu regiuni biologice scurte și comparabile; evitați matrici Python uriașe pe cromozomi întregi.

Valorile implicite diferă: demo +1/−1/−1, NW +1/−1/−2, SW +3/−3/−2. Pentru verificare față de Biopython, setați aceleași valori și documentați-le. Nu comparați direct scoruri obținute cu reguli diferite.

`pairwise2` este depreciat, dar funcționează în mediul cursului. Pentru cod nou, alternativa recomandată este `PairwiseAligner`. TODO-urile NW/SW sunt intenționate.

## Instrumente externe

[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) și [Clustal Omega](https://www.ebi.ac.uk/jdispatcher/msa/clustalo) sunt folosite în browser. `blastn`, `clustalo`, `needle` și `water` nu sunt incluse ca executabile în imaginea actuală. Instalarea lor nu este necesară pentru practica Python.

Setul proteic corectat conține p53 de om, șoarece și pește zebră; consultați [proveniența](../../data/sample/lab12-provenance.md). Nu aplicați demo-ul de distanțe pentru nucleotide unui FASTA proteic.

## Predare

[assignment.md](assignment.md) este lista completă de cerințe și punctaj. Toate fișierele sunt în `labs/02_alignment/submissions/<handle>/`, inclusiv `notes.md` și instrucțiunile de reproducere. PR-ul rămâne în repository-ul privat; în LMS trimiteți linkul și SHA-ul. Nu se cere suplimentar ZIP/PDF pentru aceeași lucrare.

Pentru actualizări ale materialelor, urmați [ghidul de actualizare](../../docs/update-course-materials.md); nu înlocuiți propriile submissions și nu combinați istoricul repository-ului public cu cel privat.
