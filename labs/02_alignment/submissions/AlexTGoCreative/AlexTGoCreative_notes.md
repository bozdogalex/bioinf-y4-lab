# Săptămâna 2 — Sequence Alignment  
Autor: AlexTGoCreative

## Date folosite
Am folosit fișierul `my_tp53.fa` din `data/work/AlexTGoCreative/lab01/`, care conține două secvențe:
- `TP53_human`
- `TP53_mouse`

Am rulat următoarele scripturi:
```bash
python ex01_global_nw.py --fasta ../../data/sample/AlexTGoCreative/lab01/my_tp53.fa --i1 0 --i2 1
python ex02_global_nw.py --fasta ../../data/sample/AlexTGoCreative/lab01/my_tp53.fa --i1 0 --i2 1


Rezultate și observații

Alinierea globală (Needleman–Wunsch) a arătat similarități pe întreaga lungime a secvențelor.

Alinierea locală (Smith–Waterman) a evidențiat doar regiunile cele mai conservate, cu un scor mai mare în zona centrală a genei.

Reflecție

Alinierea globală este preferată atunci când cele două secvențe sunt de lungimi similare și dorim o comparație completă (ex: două gene ortoloage).

Alinierea locală este preferată atunci când doar o parte a secvențelor este similară, de exemplu pentru identificarea unui domeniu proteic comun sau a unei regiuni conservate într-o genă mare.