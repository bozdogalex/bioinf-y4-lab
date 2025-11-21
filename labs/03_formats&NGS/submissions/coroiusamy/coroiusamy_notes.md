# Laborator 03: Note și Reflecție

## 1. Fișier FASTQ Utilizat

Pentru acest laborator, am folosit o rulare de secvențiere (Run) din **European Nucleotide Archive (ENA)**, relevantă pentru studii asociate cu gena **TP53**.

- **Accession ID (ENA/SRA):** [ERR179724](https://www.ebi.ac.uk/ena/browser/view/ERR179724)  
- **Link către ENA:** https://www.ebi.ac.uk/ena/browser/view/ERR179724  

### Statistici QC obținute:
- **Total Reads:** 460,577  
- **Mean Read Length:** 99.53  
- **N-Rate:** 0.0000  
- **Mean Phred Score:** 34.24  

---

## 2. Reflecție: De ce este esențială verificarea calității (QC) înainte de analiza variantelor?

Verificarea calității (**QC**) datelor FASTQ este un prim pas fundamental și **non-negociabil** în orice analiză bioinformatică, în special înainte de identificarea variantelor genetice (*variant calling*).  
Principiul de bază este **„Garbage In, Garbage Out (GIGO)”**.

### QC-ul este esențial din următoarele motive:

####  1. Prevenirea variantelor fals-pozitive
Un scor **Phred** scăzut (ex. Q10) indică o probabilitate mare (1 din 10) ca o bază să fie identificată greșit de secvențiator.  
Dacă o astfel de eroare de secvențiere (de ex., un **‘T’** citit greșit unde ar fi trebuit să fie un **‘C’**) este confundată cu o mutație biologică reală (un SNP), vom raporta o **variantă falsă**, care nu există în proba biologică.

####  2. Prevenirea variantelor fals-negative
Citirile (**reads**) de calitate foarte slabă sau cu un **N-rate** mare (multe baze ‘N’) s-ar putea să nu se alinieze (*mapeze*) corect pe genomul de referință sau ar putea fi filtrate complet de algoritmul de aliniere.  
Dacă o variantă biologică reală este acoperită doar de astfel de citiri de proastă calitate, vom rata complet detectarea acelei variante.

####  3. Eliminarea zgomotului tehnic (adaptoare)
Procesul de secvențiere poate include citirea **secvențelor de adaptor** folosite în pregătirea probelor.  
Dacă aceste secvențe tehnice (non-biologice) nu sunt identificate și eliminate (un pas numit **trimming**, decis pe baza QC-ului), ele pot cauza alinierea incorectă a citirii, generând **inserții sau deleții false** în raportul final de variante.

---

###  Concluzie
Pe scurt, **QC-ul** ne asigură că datele care intră în aliniere și variant calling sunt **de încredere** și reflectă **biologia probei**, nu erorile tehnice ale secvențiatorului.
