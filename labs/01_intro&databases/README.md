# Laboratorul 1 Baze biologice și baze de date

Începem cu [primerul de genetică și genomică](../../docs/genetics-genomics-primer-ro.md): celulă, ADN, cromozom, genă, genom, ARN și proteină. Exemplul unei gene reale apare după aceste noțiuni.

## Mediu și repository

Lucrați în repository-ul vostru **privat**, partajat cu `bozdogalex`. [Predare](../../docs/git-workflow.md) și [actualizarea unui repository existent](../../docs/update-course-materials.md). În Codespaces, deschideți Terminal → New Terminal și rulați din rădăcina repository-ului:

```bash
python labs/00_smoke/smoke.py
python -c "import Bio; print(Bio.__version__)"
```

Primul script trebuie să afișeze `ok`. Python și Biopython sunt deja în imagine. Jupyter este opțional. Puneți între ghilimele toate căile cu `&`, inclusiv în Bash.

## Demonstrații ghidate

```bash
python "labs/01_intro&databases/demo02_seq_ops.py"
python "labs/01_intro&databases/demo01_entrez_brca1.py"
```

Prima demonstrație arată transcrierea, traducerea unei secvențe artificiale, reverse complement și GC. `*` reprezintă un stop, iar pozițiile sunt indexate de la zero. A doua citește offline transcriptul BRCA1 **NM_007294.4**, de 7088 nucleotide. Nu este întregul locus genomic și nu se traduce automat de la prima bază; regiunea CDS este adnotată în GenBank.

Pentru demonstrația opțională de acces NCBI, configurați emailul real:

```bash
read -r -p "Email pentru NCBI: " NCBI_EMAIL
export NCBI_EMAIL
python "labs/01_intro&databases/demo01_entrez_brca1.py" --refresh
```

Descărcarea folosește exact aceeași accesie și salvează implicit în `data/work/demo/lab01/brca1.gb`. Dacă serviciul nu răspunde, reluați fără `--refresh` pentru copia inclusă. `demo03_dbsnp.py` este o extensie de depanare: configurați emailul din mediu înainte de rulare, conform [provocărilor pentru studenți](../../docs/lab12-student-challenges.md).

## Exercițiul de descărcare și GC

Copiați scheletul în folderul personal, păstrând originalul:

```bash
read -r -p "GitHub handle: " HANDLE
mkdir -p "labs/01_intro&databases/submissions/$HANDLE"
cp "labs/01_intro&databases/ex01_multifasta_gc.py" "labs/01_intro&databases/submissions/$HANDLE/"
```

Completați TODO-urile pentru Entrez, descărcare, citire FASTA și afișare GC. Scheletul oprește intenționat execuția până este completat. Analiza GC acceptă nucleotide, nu proteine. După implementare:

```bash
python "labs/01_intro&databases/submissions/$HANDLE/ex01_multifasta_gc.py" \
  --email "$NCBI_EMAIL" --accession NM_000546.6 \
  --out "data/work/$HANDLE/lab01/human_tp53.fa"
```

O accesie produce o singură înregistrare. Pentru comparația din Lab 2, descărcați și NM_011640.3 și NM_131327.2 în fișiere distincte, apoi reuniți-le într-un multi-FASTA. Verificați organismul și tipul moleculei pentru fiecare. În sesiunea ghidată puteți folosi `data/sample/tp53_dna_multi.fasta`; pentru evaluare documentați propriile descărcări și proveniența lor. Păstrați datele de lucru în `data/work/$HANDLE/lab01/`, ignorat de Git.

## Livrabile

În PR-ul privat includeți copia completată a exercițiului și `labs/01_intro&databases/submissions/<handle>/notes.md`, cu comenzile, accesii și versiuni, rezultate GC, convenția pentru baze ambigue și o observație despre transcript versus genă. Adăugați propria linie în roster numai în repository-ul privat. Menționați dacă demonstrația NCBI a fost live sau offline; extensia dbSNP nu blochează laboratorul de bază. Trimiteți linkul PR și SHA-ul în LMS.

Continuați cu [Lab 2](../02_alignment/README.md). Pentru început, o secvență artificială scurtă este suficientă ca să înțelegem algoritmul; datele biologice vin după aceea.
