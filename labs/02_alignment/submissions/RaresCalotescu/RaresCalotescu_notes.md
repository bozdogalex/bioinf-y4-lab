***demo01***
[INPUT]
A: CTCAAAA
B: TTTCCCC

[GLOBAL] top alignment:
CT-CAAAA
.| | ...
TTTC-CCC
  Score=-4

[LOCAL] top alignment:
2 TC
  ||
3 TC
  Score=2

***demo02***
NM_000546.6-NM_011640.3,1286,0.7221,1781
NM_000546.6-NM_131327.2,1681,0.7528,2233
NM_011640.3-NM_131327.2,1293,0.7260,1781

Alinierea globală este preferată atunci când cele două secvențe sunt similare pe întreaga lor lungime (de exemplu, variante ale aceleiași gene).
Alinierea locală este preferată când secvențele diferă mult ca lungime sau doar fragmente scurte sunt omoloage (de exemplu, domenii proteice comune între proteine diferite).


Alinierea globală este potrivită când cele două secvențe au aproximativ aceeași lungime și se dorește compararea lor complete.
Alinierea locală se folosește când doar o porțiune a secvențelor este similară, fiind utilă pentru identificarea regiunilor omoloage sau a motivelor conservate.