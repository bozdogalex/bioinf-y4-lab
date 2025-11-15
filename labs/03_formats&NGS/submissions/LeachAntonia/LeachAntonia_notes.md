1. Pentru exercițiu am folosit următorul fișier FASTQ:

Accesssion SRA: SRR292241
Link ENA: https://www.ebi.ac.uk/ena/browser/view/SRR292241

2. De ce este esențială verificarea calității datelor înainte de analiza variantelor?

Verificarea calității datelor (QC) este un pas critic înainte de analiza variantelor deoarece:

Datele de calitate slabă introduc erori în identificarea variantelor.
Citirile cu scoruri Phred mici, baze “N” sau erori de secvențiere pot genera variante false pozitive sau pot ascunde variante reale.

Citirile contaminate sau netăiate afectează alinierea.
Prezența adaptorilor, citirile prea scurte sau secvențele contaminate duc la alinieri greșite, care la rândul lor distorsionează analizele ulterioare.

QC-ul previne irosirea resurselor de calcul.
Analiza variantelor este costisitoare computațional. Dacă datele sunt proaste, pipeline-ul poate eșua sau produce rezultate inutilizabile.

QC-ul permite stabilirea celor mai bune pași de procesare.
În funcție de rezultatele QC, pot fi aplicate trimming, filtrare sau chiar re-secvențiere dacă datele sunt prea slabe.

Concluzie:
Verificarea calității este esențială pentru a garanta că datele sunt suficient de curate, corecte și fiabile înainte de a continua cu analiza variantelor, care depinde puternic de precizia fiecărei citiri.