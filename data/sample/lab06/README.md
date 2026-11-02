# Date sintetice pentru laboratoarele 6–7

`expression_matrix.csv`: 120 de gene fictive × 48 de probe sintetice independente. Valorile sunt expresie continuă pregătită pe o scară asemănătoare log2. Nu sunt counts, probe de pacienți sau rezultate experimentale.

`sample_metadata.csv`: 24 de controale și 24 de probe tratate simulate. Aceste etichete nu sunt folosite la calcularea corelațiilor sau la Louvain.

`generation_groups.csv`: grupurile generatorului, folosite doar pentru a compara rezultatul inferenței cu mecanismul simulării. P1 și P3 urmăresc variabile latente opuse; P2 urmărește o variabilă latentă independentă; 30 de gene au numai zgomot.

`provenance.json`: parametrii și avertismentele simulării. Regenerare din rădăcina repository-ului:

```bash
python labs/06_wgcna/generate_teaching_data.py
```

Seed-ul și versiunile bibliotecilor trebuie păstrate pentru reproducere. Fișierele CSV livrate reprezintă intrarea stabilă a laboratorului; nu este nevoie să le regenerați înainte de rulare. Identificatorii G001–G120 nu se trimit la servicii de enrichment sau de adnotare biologică.

Date și generator pregătite cu asistență AI, sub îndrumarea cadrului didactic.
