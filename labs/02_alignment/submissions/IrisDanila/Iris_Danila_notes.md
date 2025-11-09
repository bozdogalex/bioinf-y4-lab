# Săptămâna 2 — Sequence Alignment  
Autor: IrisDanila

## Date folosite
Am folosit fișierul `tp53_mrna.fa` din `data/work/IrisDanila/lab01/`, care conține cinci secvențe TP53 mRNA:
- `NM_001126117.2` (TP53 transcript variant 7, 2063 bp)
- `NM_001126115.2` (TP53 transcript variant 5, 2003 bp)

Am rulat următoarele scripturi:
```bash
python ex01_global_nw.py --fasta data/work/IrisDanila/lab01/tp53_mrna.fa --i1 0 --i2 1
python ex02_local_sw.py --fasta data/work/IrisDanila/lab01/tp53_mrna.fa --i1 0 --i2 1

Rezultate și observații

Alinierea globală (Needleman–Wunsch) a arătat similarități pe întreaga lungime a secvențelor. Primele 500 bp sunt identice (scor: 500), indicând o regiune 5' UTR conservată între variante.

Alinierea locală (Smith–Waterman) a evidențiat regiunile cele mai conservate, obținând un scor de 5889. Am observat o regiune cu gap-uri (poziția ~650-680), sugerând diferențe de splicing alternativ între cele două variante de transcript.

Reflecție
Alinierea globală este preferată atunci când cele două secvențe sunt de lungimi similare și dorim o comparație completă (ex: compararea variantelor de transcript ale aceleiași gene, detectarea SNP-urilor).

Alinierea locală este preferată atunci când doar o parte a secvențelor este similară, de exemplu pentru identificarea domeniilor proteice conservate, căutarea de motive funcționale, sau când lungimile secvențelor diferă semnificativ.