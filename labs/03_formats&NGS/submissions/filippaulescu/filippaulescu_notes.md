Laborator 03: QC și Reflecție

1. Fișierul Meu

Am folosit date din ENA pentru a studia gena TP53.

Accession: ERR179724

Fișier descărcat: your_reads.fastq.gz (5.1 milioane de citiri)

Statistici QC (Rezultat ex04_fastq_qc.py):

Datele obținute arată o calitate excelentă:

Total Reads: 5,092,543 (Volum mare de date).

Lungime Medie: 75.00 (Citiri scurte, tipice Illumina).

Rata N: 0.0000 (Perfect. Fără baze necunoscute).

Scor Phred Mediu: 34.47 (Calitate excelentă, P > 99.96%).

2. Reflecție: De ce e important QC-ul?

QC (Quality Control) este primul pas și este OBLIGATORIU. Fără el, analiza ulterioară este inutilă. Regula e simplă: Garbage In, Garbage Out (GIGO).

De ce facem QC (motivul principal)

QC-ul ne asigură că variațiile pe care le găsim (mutațiile) sunt biologice, nu erori de mașină (artefacte).

Problema: Scor Phred Scăzut

Consecința: Găsim variante fals-pozitive (o eroare de citire e confundată cu o mutație reală).

Soluția QC: Ne permite să tăiem (trim) bazele de calitate slabă.

Problema: Baze 'N' sau Adaptoare

Consecința: Mapare (aliniere) greșită sau ratarea mutațiilor (fals-negative).

Soluția QC: Curățăm adaptorii și reads-urile neclare pentru o aliniere precisă.

Concluzie

Datele mele (Phred 34.47, N-Rate 0.0000) sunt curate. QC-ul confirmă că pot trece la aliniere și variant calling cu încredere.