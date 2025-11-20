FastQ folosit:

-Accession: ERR071289  
-Sursă: ENA (European Nucleotide Archive)  
-Cale locală: `data/work/rauldavid35/lab03/ERR071289.fastq.gz`  
-Link (info): https://www.ebi.ac.uk/ena/browser/view/ERR071289



De ce este esențială verificarea calității datelor înainte de analiza variantelor?

Pe scurt: dacă pornesc cu citiri „murdare”, tot ce urmează (mapare, VCF) se strică.  
QC-ul îmi spune din start dacă datele sunt ok sau trebuie să tai/adaptez.
De exemplu:
-Phred mic → baze nesigure, apar “variante” false.
-N-uri/adaptoare → mapare mai slabă, aliniamente ambigue.
-Acoperire neuniformă (bias GC) → variante reale pot fi ratate.
-Duplicate/artefacte → acoperire umflată, filtre înșelate.