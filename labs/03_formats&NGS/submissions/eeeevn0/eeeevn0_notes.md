
**Daniela Ispas**
**Eugenia Vacarciuc**


**Fișierul FASTQ utilizat**

Pentru acest laborator am descărcat fișierul FASTQ cu accession **SRR000001** din **ENA (European Nucleotide Archive)**, folosind scriptul `ex01_fetch_fastq.py`. 

Link:[SRR000001.fastq.gz — ENA FTP](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/001/SRR000001/SRR000001.fastq.gz)  

Fișierul conține citiri brute de ADN (reads) generate prin tehnologia **Illumina** și a fost salvat local în `data/work/eeeevn0/lab03/SRR000001.fastq.gz`. Ulterior, același fișier a fost analizat pentru verificarea calității (QC).

**Reflecție — de ce este importantă verificarea calității (QC) înainte de analiza variantelor?**

Verificarea calității este un pas esențial în analiza datelor de secvențiere,deoarece permite identificarea citirilor eronate, a bazelor necunoscute și a zonelor cu calitate scăzută.  
Am observat că fișierul nostru FASTQ inițial era corupt, ceea ce ar fi putut duce la erori mari în analiză. QC este esențială pentru a verifica dacă citirile sunt corecte și de bună calitate.Astfel evităm rezultate false și ne asigurăm că analiza variantelor se bazează pe date curate și de încredere. În final, QC contribuie la obținerea unor rezultate precise și valide din punct de vedere biologic.

**Rezultatele obtinute**
Reads: 236003
Mean length: 254.62
N rate: 0.0005
Mean Phred: 25.96