**Daniela Ispas**
**Eugenia Vacarciuc**

**Secvențe FASTA folosite**

Pentru acest laborator am refolosit fișierul **SRR000001.fastq.gz** de la laboratorul 3, descărcat din baza de date publică **ENA (European Nucleotide Archive)**.  
Fișierul conține citiri de ADN (reads) obținute prin tehnologia **Illumina**.  

Apoi l-am convertit în format **FASTA** pentru a putea fi folosit în analiza filogenetică și am extras trei citiri, pe care le-am salvat într-un fișier multi-FASTA `data/work/eeeevn0/lab04/my_sequences.fasta`

Cele trei secvențe FASTA folosite:

>SRR000001.1
ATTCTCCTAGCCTACATCCGTACGAGTTAGCGTGGGATTA...
>SRR000001.4
AAAGGAGTGTGACATTCTGTGTTCCACATGCATCGACTAG...
>SRR000001.5
CTCTTTCCCTCTTCCTCCTCTTCCTCCTTCTTATTCCTCT...

**Rezultate**
Număr secvențe: 3
Lungimi inițiale: [266, 240, 276] → trunchiate la 240 baze pentru a fi comparabile

Matricea de distanțe (identity):
SRR000001.1  0.000000
SRR000001.4  0.775000    0.000000
SRR000001.5  0.700000    0.758333    0.000000

Arbore filogenetic (Newick): 
(SRR000001.1:0.35833,SRR000001.4:0.41667,SRR000001.5:0.34167)Inner1:0.00000;
  

**Reflecție**

Matricea de distanțe ne arată cât de diferite sunt secvențele între ele,
dar arborele filogenetic merge mai departe, el evidențiază relațiile de înrudire și modul în care aceste secvențe s-au ramificat dintr-un strămoș comun.

În cazul nostru, SRR000001.1 și SRR000001.5 sunt mai apropiate între ele decât de SRR000001.4,
ceea ce sugerează existența unui strămoș comun mai recent.

Astfel, nu observăm doar diferențele numerice dintre secvențe, ci și legăturile evolutive care le unesc. Analiza filogenetică transformă datele brute într-o imagine clară și intuitivă a relațiilor dintre secvențe.
