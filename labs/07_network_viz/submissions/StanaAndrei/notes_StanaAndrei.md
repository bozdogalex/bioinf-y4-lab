Iată răspunsul pe scurt:

### 1. Metoda de layout folosită
Am folosit **Spring Layout** (algoritmul Fruchterman-Reingold).
* **Cum funcționează:** Simulează un sistem fizic în care nodurile se resping (ca magneții), iar muchiile le atrag (ca niște arcuri).
* **Efect:** Genele strâns conectate se atrag automat, formând clusteri vizuali distincți care corespund modulelor.

### 2. Avantaje vizualizare (Lab 7) vs. Analiză numerică (Lab 6)
* **Validare intuitivă:** Poți verifica imediat dacă matematica a funcționat: nodurile cu aceeași culoare (același modul) ar trebui să fie grupate spațial, nu împrăștiate.
* **Identificarea rolurilor:** Distingi între **hub-uri centrale** (esențiale pentru un modul) și **hub-uri "punte"** (care fac legătura între două module diferite), lucru imposibil de văzut într-un tabel Excel.
* **Detectarea zgomotului:** Vezi dintr-o privire dacă rețeaua este un "hairball" (prea densă/zgomotoasă) sau dacă are o structură biologică clară.

