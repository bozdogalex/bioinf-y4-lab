# Saptamana 2 — Sequence Alignment Notes (Papaia462819)

- Date folosite: `my_tp53.fa` (secventele 0 si 1) si `nm000546.fa` pentru verificari suplimentare.
- Scoruri utilizate in implementari: match=1, mismatch=-1, gap=-2 pentru NW; match=3, mismatch=-3, gap=-2 pentru SW.
- Implementarile sunt rulate pe subseturi din acele fisiere pentru a testa atat aliniere globala, cat si locala.

## Reflectie
Alinierea globala (Needleman-Wunsch) este preferata cand secventele au lungimi comparabile si ne asteptam ca intreaga regiune sa fie omologa (ex. doua izoforme apropiate). Alinierea locala (Smith-Waterman) este mai potrivita cand doar un segment scurt este similar sau cand secventele difera mult ca lungime si vrem sa identificam domenii conservate.
