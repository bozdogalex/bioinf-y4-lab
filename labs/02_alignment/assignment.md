# Laboratorul 2 Cerințe și evaluare

Lucrați în propriul repository privat, partajat cu `bozdogalex`. Această pagină este lista de evaluare pentru Lab 2; [README](README.md) descrie demonstrațiile. Predați în LMS **URL-ul PR-ului privat și SHA-ul complet**. Nu se cere în paralel o arhivă ZIP sau un raport PDF.

## Date și reproducere

Folosiți cel puțin trei secvențe comparabile descărcate în Lab 1: de exemplu transcripturi ale aceleiași gene din organisme diferite. Specificați accesii cu versiuni, organisme, tipul moleculei, data și sursa descărcării. Un transcript, o secvență genomică și o proteină nu sunt interschimbabile.

Fișierele comune și perechea artificială sunt permise pentru demonstrații și depanare. Pentru evaluare folosiți și documentați datele proprii. Păstrați datele mari în `data/work/<handle>/lab01/` și includeți instrucțiuni de descărcare reproductibilă. Nu predați date personale sau sensibile.

Puneți codul, `README.md` cu comenzile și versiunile, `notes.md` (aproximativ o pagină), matricea de distanțe și extrasul MSA în `labs/02_alignment/submissions/<handle>/`. Pentru comparații foarte mari, salvați un extras și descrieți cum se regenerează rezultatul complet.

## Task 1 Distanțe pe poziții aliniate 3 puncte

Calculați p-distance pentru fiecare pereche din set. Aliniați secvențele înainte de numărare sau folosiți un MSA. Definiți tratamentul gap-urilor și bazelor ambigue și afișați numărul pozițiilor comparate; fără poziții comparabile, raportați `NA`, nu zero. Hamming este adecvat doar pentru șiruri de aceeași lungime pe poziții deja corespunzătoare. Trunchierea nu înlocuiește alinierea.

Produceți o matrice sau tabelul perechilor. În `notes.md`, identificați perechea cu distanța cea mai mică și discutați limita interpretării. O distanță observată pe un fragment nu demonstrează singură o relație filogenetică.

## Task 2 Aliniere globală și locală 4 puncte

Completați TODO-urile din copiile `ex01_global_nw.py` și `ex02_local_sw.py`. Verificați întâi pe perechi artificiale scurte, apoi pe două regiuni biologice comparabile, de aproximativ 50–200 nucleotide. Înregistrați accesii, coordonate și motivul alegerii regiunii. Pentru coordonate Python specificați convenția zero-based, cu capătul drept exclus.

Rulați și Biopython pe aceleași secvențe, cu aceleași scoruri pentru match, mismatch și gap. Comparați scorul și validitatea alinierii; pot exista trasee optime diferite. Nu folosiți `globalxx`/`localxx` cu gap-uri gratuite ca unic exemplu biologic.

Explicați diferența global/local, lungimile regiunilor aliniate și efectul gap-urilor. Includeți un fragment ilustrativ. Alocare: 2p implementare și verificare, 2p comparație și interpretare.

## Task 3 Aliniere multiplă 3 puncte

Aliniați cele trei secvențe cu Clustal Omega sau un instrument echivalent. Alegeți tipul corect de moleculă. Salvați rezultatul text și includeți un extras relevant, nu doar un link temporar al serviciului.

Marcați o regiune conservată, separați observația de ipoteza funcțională și explicați ce aduce a treia secvență față de comparațiile pe perechi. Dacă folosiți proteine, toate intrările trebuie să reprezinte familia de proteine aleasă, cu organismele și isoformele verificate.

## Bonus Aliniere semi globală 1 punct

Configurați o aliniere semi-globală, precizați care capete sunt nepenalizate și motivați un caz de utilizare. Bonusul nu înlocuiește cerințele de bază.

## Integritate și colaborare

Puteți lucra individual sau în perechi conform [politicii cursului](../../docs/policies.md). Identificați ambii autori în raport și PR dacă lucrați în pereche. Menționați sursele externe și asistența AI în README. Testele de mediu verzi nu reprezintă evaluarea corectitudinii științifice.
