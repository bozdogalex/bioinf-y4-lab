# Săptămâna 2 — Assignment (Sequence Alignment)

## Instrucțiuni generale
- Folosiți **DOAR** propriile secvențe descărcate în Lab 1 din NCBI (stocate local în `../../data/work/<handle>/lab01/`).
- Puteți lucra în Jupyter (notebook) sau fișiere `.py`.
- **Predare (Moodle)**: încărcați un fișier ZIP numit `lab02_alignment_<handle>.zip` care conține:
  - scripturile / notebook-urile voastre;
  - un fișier `README.txt` (max 10 rânduri) cu **pașii de rulare**, versiunea de Python și dependențele;
  - un fișier `notes.pdf` (max 1 pagină) cu răspunsurile / interpretările cerute mai jos.

---

## Task 1 — Distanțe perechi (3p)
- Implementați și rulați **Hamming** (numai pentru perechi de **aceeași lungime**) **sau** **p-distance** (proporția pozițiilor diferite) pentru toate perechile dintr-un subset de **≥3 secvențe** din fișierul vostru.
- Produceți o **matrice de distanțe** (triunghiul superior este suficient).
- În `notes.pdf` (2–3 rânduri): **care două secvențe sunt cele mai apropiate** și **de ce** (argument biologic / logic).
*Hint:* dacă lungimile diferă, **nu** folosiți Hamming. Alegeți: (a) p-distance pe aliniamente pairwise brute (`globalxx` în Biopython) sau (b) trunchiați la lungimea minimă și **motivați alegerea**.

---

## Task 2 — Pairwise alignments (4p)
- Alegeți două secvențe din dataset (specificați ID-urile).
- Rulați două aliniamente pairwise cu Biopython:
  - global (ex. `pairwise2.align.globalxx` sau variantă cu scor match/mismatch/gap),
  - local (ex. `pairwise2.align.localxx`).
- În `notes.pdf` (max 6–7 rânduri):
  - Comparați **global vs. local** (regiuni aliniate, număr/poziție gap-uri, scor).
  - Includeți un **mic fragment** din aliniere unde local găsește o potrivire pe care global o “forțează” cu gap-uri (dacă identificați un astfel de caz).

---

## Task 3 — MSA online (3p)
- Alegeți ≥3 secvențe din dataset.
- Rulați aliniere multiplă (MSA) cu **Clustal Omega (EBI)** sau un alt instrument online echivalent.
- Exportați rezultatul MSA (text) și includeți un **extras relevant** în `notes.pdf` (sau link permanent, dacă există).
- În `notes.pdf` (max 6–7 rânduri):
  - Marcați o **regiune conservată** (motif / segment identic) și explicați de ce credeți că este conservată.
  - Comentați când **MSA ajută** interpretarea comparativ cu aliniamentele pe perechi.

---

## Bonus — Semiglobal (+1p)
- Rulați un **aliniament semiglobal** (implementare simplificată sau setare dintr-un tool care nu penalizează gap-urile la capete).
- În `notes.pdf` (3–4 rânduri): **când** ați prefera semiglobal vs. global/local?

---

## Punctaj & criterii
- Task 1: **3p** — corectitudinea calculelor + matrice clară.
- Task 2: **4p** — rulare corectă + **interpretare** global vs. local cu fragment exemplu.
- Task 3: **3p** — MSA corectă + **identificare motiv conservat** + comparație cu pairwise.
- Bonus: **+1p** — scenariu semiglobal motivat.

---

## Integritate academică 
- Toate lucrarile pot fi executate singur, sau in perechi, conform [docs/policies.md](../../docs/policies.md)  
- Dacă lucrați în perechi, indicați ambele nume în `notes.pdf` și includeți ambii autori în arhiva `.zip`.  
- Dacă folosiți resurse externe (inclusiv AI), notați sursa pe scurt în `README.txt`.