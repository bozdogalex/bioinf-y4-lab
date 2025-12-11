# Laborator 4 — Filogenetică
---

## 🔹 Setul de date utilizat
Am folosit o rulare de secvențiere din **European Nucleotide Archive (ENA)**, asociată cu studii asupra genei **TP53**:

- **Run ID (ENA/SRA):** [ERR179724](https://www.ebi.ac.uk/ena/browser/view/ERR179724)  
- **Descriere:** date de tip *Illumina paired-end sequencing reads*, extrase și convertite din FASTQ în FASTA.  
- **Fișierul multi-FASTA folosit:** `data/work/coroiusamy/lab04/my_sequences.fasta`  
- **Conținut (3 secvențe):**
  ```
  >ERR179724.1
  ATACGTAGGTGGCGAGCGTTGTCCGGATTTACTGGGCGTAAAGGGAGCGTAGGCGGACTT...
  
  >ERR179724.2
  ATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCGCGCAGGTGGTTAA...
  
  >ERR179724.3
  ATACGTAGGTGGCAAGCGTTATCCGGAATCATTGGGCGTAAAGGGTGCGTAGGTGGCGTA...
  ```

---

## 🔹 Rezultate
- **Matricea de distanțe (p-distance / identity):**
  ```
  ERR179724.1  0.000000
  ERR179724.2  0.274510    0.000000
  ERR179724.3  0.294118    0.215686    0.000000
  ```
- **Arbore Neighbor-Joining (Newick):**
  ```
  (ERR179724.1:0.17647,ERR179724.2:0.09804,ERR179724.3:0.11765)Inner1:0.00000;
  ```
- **Vizualizare ASCII:**
  ```
    _________________________________________________________________ ERR179724.1
   |
   |___________________________________ ERR179724.2
   |
   |__________________________________________ ERR179724.3
  ```

---

## 🔹 Reflecție — Ce informații suplimentare oferă arborele filogenetic față de o simplă matrice de distanțe?

O **matrice de distanțe** ne arată doar *cât de diferite* sunt două secvențe între ele, dar nu oferă o imagine clară despre **relațiile evolutive** dintre toate secvențele simultan.  
În schimb, un **arbore filogenetic**:

- vizualizează *ierarhia* similarităților dintre secvențe;  
- evidențiază care secvențe au un **strămoș comun mai apropiat**;  
- arată posibile **grupuri (clustere)** de secvențe înrudite;  
- permite inferențe evolutive și corelarea cu divergența funcțională sau geografică.

Astfel, arborele nu doar cuantifică diferențele, ci le **structurează evolutiv**, oferind o perspectivă mai interpretabilă asupra relațiilor dintre secvențe.

