# Introducere în genetică și genomică pentru bioinformatică

Acest material este un suport de predare pentru studenți care cunosc programare, dar au puține cunoștințe de biologie. Poate fi parcurs înaintea alinierii secvențelor, ca extindere a primelor cinci slide-uri. Explicațiile sunt formulate pentru prezentare orală; întrebările scurte sunt pentru clasă. Exemplele de secvență fără accesie sunt artificiale.

Firul lecției este: **ce există în celulă, cum este organizată informația genetică, cum este folosită și cum ajunge într-un fișier pe care îl analizăm**. O genă concretă, precum BRCA1 sau TP53, apare abia la aplicarea noțiunilor. Nu este nevoie ca studenții să învețe mai întâi funcția unei gene pentru a înțelege vocabularul.

## 1 De ce începem cu biologia

În laborator vom vedea multe șiruri de litere. Le putem citi din fișiere, le putem compara și putem calcula statistici. Dar înainte de asta trebuie să înțelegem ce reprezintă. Același calcul poate avea sens pentru un tip de date și poate fi lipsit de sens pentru altul.

Să pornim de la o întrebare generală. O celulă musculară și un neuron din corpul aceleiași persoane au funcții și forme foarte diferite. De unde vine această diferență? Au nevoie de un set complet diferit de gene?

În general, celulele somatice ale unei persoane au aproape același genom nuclear. Diferențele dintre tipurile de celule depind în mare măsură de genele exprimate, de cantitatea produselor lor și de felul în care activitatea lor este reglată. Există și excepții biologice, dar nu avem nevoie de ele pentru a înțelege primul principiu: informația disponibilă și informația folosită într-un anumit moment sunt lucruri distincte.

Pentru un informatician, este tentant să compare genomul cu un program. Analogia poate ajuta la ideea de informație, dar are limite. Celula este un sistem fizic și chimic, în care moleculele interacționează, iar funcționarea depinde de mediu, structură și cantitate. Nu există un singur procesor care execută genomul linie cu linie.

Astăzi vom construi vocabularul necesar pentru a trece de la moleculă la reprezentarea ei în date. La fiecare pas, încercăm să răspundem la două întrebări: ce este obiectul biologic și cum îl reprezentăm în calculator?

*Întrebare inițială: dacă două celule au aproape același ADN, ce alte date am putea măsura ca să înțelegem de ce funcționează diferit? Acceptă răspunsuri intuitive despre activitate, ARN și proteine; revino la ele la secțiunea despre expresie.*

## 2 Celula și locul materialului genetic

Celula este unitatea de bază a organismelor vii. Are o membrană care delimitează un interior și conține molecule și structuri care permit transformarea energiei, sinteza componentelor și răspunsul la mediu.

La organismele eucariote, precum animalele, plantele și fungii, cea mai mare parte a ADN-ului se află în nucleu. Există și ADN în mitocondrii, iar la plante și în cloroplaste. Prin urmare, când spunem „ADN-ul unei celule”, trebuie să fim atenți dacă discutăm genomul nuclear sau și genomurile organitelor.

Bacteriile nu au un nucleu delimitat prin membrană. Materialul lor genetic este organizat diferit și pot avea și plasmide, molecule de ADN distincte de cromozomul principal. Nu trebuie să transferăm automat schema celulei umane la toate organismele.

Pentru laborator, este suficient să rețineți relația dintre niveluri: organismul este alcătuit din celule, celulele conțin material genetic, iar materialul genetic are o organizare și poate fi măsurat experimental.

Mai există o distincție: o probă biologică nu conține neapărat un singur tip de celulă sau un singur organism. O probă de țesut poate include tipuri celulare diferite. O probă de mediu poate conține microorganisme diferite. Când ajungem la fișier, aceste informații despre probă nu dispar din importanță.

Aceasta este una dintre explicațiile pentru care bioinformatica folosește atât secvențe, cât și metadate. Literele singure nu ne spun dacă provin din sânge, dintr-o cultură bacteriană sau dintr-un amestec de organisme.

*Verificare: un fișier cu ADN uman provine obligatoriu doar din nucleu? Răspunsul trebuie să amintească ADN-ul mitocondrial, fără a transforma discuția într-un inventar de excepții.*

## 3 Din ce este alcătuit ADN-ul

ADN înseamnă acid dezoxiribonucleic. Este o moleculă alcătuită din unități numite nucleotide. Fiecare nucleotid conține un zahăr, o grupare fosfat și o bază azotată. În reprezentarea unei secvențe, folosim litera bazei: A pentru adenină, C pentru citozină, G pentru guanină și T pentru timină.

Așadar, când vedeți `ACGT`, nu aveți patru gene. Aveți reprezentarea a patru nucleotide într-o anumită ordine. Informația pe care o analizăm în aceste laboratoare este, în primul rând, această ordine.

ADN-ul celular este de obicei bicatenar: are două catene. Bazele se împerechează complementar, A cu T și C cu G. Structura este cunoscută ca dublu helix. Scheletul zahăr–fosfat oferă continuitatea catenei, iar bazele sunt partea variabilă pe care o codificăm prin litere. [NHGRI despre ADN](https://www.genome.gov/about-genomics/fact-sheets/Deoxyribonucleic-Acid-Fact-Sheet)

Termenul „pereche de baze”, prescurtat bp, se referă la două baze împerecheate din catene opuse. Pentru o moleculă bicatenară de 100 bp, fiecare catenă are 100 nucleotide. În fișiere reprezentăm frecvent o singură catenă în orientarea convențională, deoarece cealaltă poate fi dedusă prin complementaritate.

O secvență este ordonată. `AAGC` și `AGAC` conțin aceleași tipuri și numere de baze, dar în altă ordine. Fracția GC este 0,5 pentru ambele, însă secvențele nu sunt identice. Această diferență dintre compoziție și ordine va motiva alinierea.

În datele reale pot apărea și simboluri ambigue. De exemplu, N indică o bază nespecificată. N nu este o a cincea bază standard. Un program trebuie să decidă cum tratează incertitudinea când calculează o statistică.

*Pe tablă: scrie AAGC și AGAC. Cere numărul de G/C, apoi întreabă dacă această statistică identifică în mod unic secvența.*

## 4 Orientarea unei secvențe

Catenele de ADN au direcție. Capetele se numesc 5′ și 3′, după pozițiile atomilor din zahăr. Pentru laborator ne interesează consecința: prin convenție, secvențele sunt de regulă scrise de la 5′ la 3′.

Cele două catene ale ADN-ului sunt antiparalele. Dacă una este scrisă de la 5′ la 3′, catena complementară așezată sub ea merge de la 3′ la 5′:

```text
5′ A T G C 3′
3′ T A C G 5′
```

Dacă vrem să raportăm și a doua catenă în convenția 5′ spre 3′, trebuie să o citim în sens invers. Rezultatul este `GCAT`, numit reverse complement al lui `ATGC`.

Operația are două componente: complementaritate și inversarea ordinii. Simplul reverse al lui ATGC este CGTA; acesta nu este reverse complement. Această diferență contează când o genă este adnotată pe catena opusă sau când potrivim o citire pe o referință.

Din perspectiva programării, orientarea este o parte a semnificației datelor. Două șiruri pot reprezenta cele două orientări ale aceleiași regiuni de ADN. Dacă ignorăm această posibilitate, putem concluziona greșit că nu există asemănare.

Nu trebuie să memorați chimia detaliată a capetelor astăzi. Trebuie să știți convenția de scriere și să puteți explica rezultatul funcției `reverse_complement()` din Biopython.

*Exercițiu de 30 de secunde: pentru AAGC, complementul așezat antiparalel este TTCG, iar reverse complement în orientarea 5′–3′ este GCTT. Cere studenților să marcheze capetele înainte să răspundă.*

## 5 Cromozom genă și genom

Un cromozom este o structură în care o moleculă lungă de ADN este asociată cu proteine. La eucariote, ADN-ul este împachetat împreună cu proteine precum histonele. Forma foarte condensată, desenată frecvent ca un X, este asociată cu anumite etape ale diviziunii. Cromozomii nu stau permanent în acea formă.

O genă este o regiune de ADN a cărei informație contribuie la producerea unui ARN funcțional sau, pentru genele codante, a unei proteine prin intermediul ARN-ului. Un cromozom conține multe gene și multe alte regiuni. Gena nu este o moleculă separată lipită de cromozom: este o regiune a ADN-ului respectiv.

Genomul reprezintă întregul material genetic al organismului. Nu este doar lista genelor care codifică proteine. Include regiuni de reglare, secvențe repetitive, introni și alte regiuni. Într-o analiză trebuie să precizăm dacă discutăm genomul nuclear, mitocondrial sau ambele.

Într-o celulă somatică umană tipică, nucleul conține 46 de cromozomi organizați în 23 de perechi. Gameții au în mod obișnuit 23. Există excepții și variații, iar unele celule, precum eritrocitele umane mature, nu au nucleu. [NHGRI despre cromozomi](https://www.genome.gov/about-genomics/fact-sheets/Chromosomes-Fact-Sheet)

Doi cromozomi omologi dintr-o pereche poartă în general aceleași tipuri de gene la poziții corespunzătoare, dar secvențele lor nu sunt obligatoriu identice. Distingeți omologii de cromatidele surori, produse prin replicarea unui cromozom.

O analogie limitată: dacă genomul este o colecție de volume, cromozomii sunt volumele, iar genele sunt anumite regiuni cu informație funcțională. Dar există și informație între aceste regiuni, iar regulile de utilizare nu urmează o lectură simplă de la prima la ultima pagină.

*Verificare: dacă un fișier conține secvența unei gene, conține automat întregul cromozom? Dacă avem toate secvențele codante, avem automat întregul genom?*

## 6 Alele și moștenire

Am spus că, pentru majoritatea regiunilor autosomale, o persoană are două copii corespunzătoare: una moștenită de la fiecare părinte. Copiile pot diferi prin secvență. Versiunile alternative la un anumit locus se numesc alele. Locus înseamnă poziția sau regiunea genomică despre care vorbim.

Dacă ambele copii au aceeași variantă la locusul analizat, spunem că persoana este homozigotă pentru acea variantă. Dacă sunt diferite, spunem heterozigotă. Termenii trebuie raportați la un locus; o persoană nu este pur și simplu „heterozigotă” în toate privințele. [NHGRI despre alele](https://www.genome.gov/genetics-glossary/Allele)

Să luăm un exemplu artificial: la o poziție de pe un autosom, o copie are A și cealaltă are G. Putem reprezenta genotipul A/G. Genotipul descrie configurația genetică, în timp ce fenotipul descrie caracteristici observabile sau măsurabile. Relația dintre ele poate depinde de multe gene, de reglare și de mediu.

Nu presupunem că orice caracteristică este determinată de o singură genă și că o literă schimbată are întotdeauna un efect vizibil. „Dominant” nu înseamnă mai frecvent în populație sau biologic mai bun; se referă la modul de exprimare a unui efect într-un anumit context genetic.

În reproducerea sexuată, meioza și recombinarea contribuie la transmiterea și amestecarea materialului genetic. Nu vom construi astăzi modele de segregare, dar această idee explică de ce indivizii aceleiași specii au multe asemănări și totuși nu au secvențe complet identice.

Pentru fișierele pe care le vom analiza, consecința este că o secvență de referință nu reprezintă toate alelele tuturor persoanelor. Ea oferă un reper pentru comparație.

*Întrebare: o diferență față de referință este automat o eroare sau o cauză de boală? Cere cel puțin două alternative: variație biologică reală și eroare de măsurare.*

## 7 Replicarea și expresia sunt procese diferite

Înainte de diviziune, celula trebuie să copieze ADN-ul. Procesul se numește replicare. Complementaritatea bazelor permite folosirea fiecărei catene ca matriță pentru sinteza unei catene noi. Replicarea transmite informația genetică către celulele rezultate.

Expresia unei gene înseamnă folosirea informației ei pentru a produce un ARN și, în cazul genelor codante, o proteină. Nu este același proces cu replicarea. O celulă poate exprima gene fără să își copieze întregul genom pentru diviziune.

Revenim la neuron și celula musculară. Ele folosesc diferit informația disponibilă: unele gene sunt exprimate mai mult, altele mai puțin, iar această activitate se schimbă în timp și în funcție de semnale. Expresia nu este doar o listă fixă de gene pornite sau oprite; cantitatea produselor contează. [NHGRI despre expresia genică](https://www.genome.gov/genetics-glossary/Gene-Expression)

Reglarea implică interacțiuni între ADN, ARN, proteine și organizarea cromatinei. Unele modificări ale cromatinei influențează accesibilitatea ADN-ului fără să schimbe ordinea bazelor; acestea intră în discuția despre reglare epigenetică. Pentru moment, rețineți distincția dintre schimbarea secvenței și schimbarea modului în care este utilizată.

Într-o analiză, secvențierea ADN-ului și măsurarea ARN-urilor răspund unor întrebări diferite. Prima poate descrie variante de secvență. A doua ne poate spune ce transcripturi sunt prezente și în ce cantități, în condițiile măsurate.

Astfel, putem avea același identificator de genă în două tabele foarte diferite: unul despre secvență și altul despre expresie. Nu interpretăm coloanele doar după numele genei; verificăm ce mărime a fost măsurată.

## 8 Transcripția și tipurile de ARN

ARN înseamnă acid ribonucleic. Nucleotidele sale conțin riboză, iar alfabetul de bază folosește U, uracil, în loc de T. ARN-ul este adesea monocatenar, dar se poate plia și poate forma structuri prin împerecheri interne.

În transcripție, o regiune de ADN este folosită ca matriță pentru sinteza ARN-ului. ARN-polimeraza produce ARN în direcția 5′ spre 3′, citind matrița în direcția opusă.

Pentru o genă codantă, putem distinge catena matriță și catena codantă. Secvența ARN-ului corespunde catenei codante, înlocuind T cu U, pentru regiunea transcrisă și înainte de a discuta procesarea ei. Enzima folosește drept matriță cealaltă catenă.

Nu toate ARN-urile au rol de mesager. ARN-ul mesager, mRNA, poate fi tradus. ARN-ul ribozomal, rRNA, este o componentă structurală și catalitică a ribozomului. ARN-ul de transfer, tRNA, contribuie la asocierea codonilor cu aminoacizii. Există și ARN-uri cu rol de reglare.

Prin urmare, schema ADN → ARN → proteină este foarte utilă pentru genele codante, dar nu înseamnă că orice ARN se transformă într-o proteină. ARN-ul are propriile funcții. ARN → ADN este posibil prin transcripție inversă; schema introductivă nu enumeră toate procesele de transfer al informației.

În demonstrația Biopython, `transcribe()` aplicată unei secvențe reprezentate ca ADN codant înlocuiește T cu U. Ea nu găsește singură promotorul, exoanele sau începutul unei gene. Este o operație pe reprezentarea pe care noi i-o furnizăm.

*Întrebare: dacă un program a înlocuit toate literele T cu U, a identificat automat o genă? Explicați ce informație biologică lipsește.*

## 9 Exoni introni și transcripturi

La multe gene eucariote, informația care va apărea în ARN-ul matur este separată în regiuni numite exoni, între care se află introni. ARN-ul inițial poate conține ambele tipuri de regiuni. Prin splicing, intronii sunt eliminați, iar exonii sunt uniți.

Un exon este o regiune păstrată în transcriptul matur; nu este obligatoriu ca întreaga sa secvență să fie tradusă în proteină. Transcriptul poate conține regiuni netraduse la capetele 5′ și 3′, numite UTR. Regiunea codantă propriu-zisă este CDS.

Această distincție explică o eroare practică frecventă: descărcăm un mRNA și îl traducem din prima literă, presupunând că am obținut proteina. Dacă primele baze aparțin UTR-ului, traducerea nu pornește din locul biologic corect.

Să desenăm o genă simplificată:

```text
ADN genomic:       exon 1 — intron — exon 2 — intron — exon 3
ARN după splicing: exon 1            exon 2            exon 3
Transcript matur:  5′ UTR | regiune codantă CDS | 3′ UTR
```

Prin splicing alternativ se pot obține combinații diferite de exoni și, împreună cu alte mecanisme, transcripturi și produse proteice diferite. Aceste forme sunt numite frecvent izoforme, cu precizarea dacă vorbim despre transcript sau proteină.

Nu confundați izoforma cu alela. Alela descrie o versiune a secvenței la un locus. Izoformele pot rezulta din procesarea sau utilizarea diferită a informației aceleiași gene. Nici expresia „o genă, o proteină” nu este o regulă universală.

Pentru analiza noastră, identificatorul unui transcript și identificatorul unei proteine trebuie păstrate separat. Numele genei singur poate fi insuficient pentru reproducerea unei comparații.

*Deschide ulterior recordul GenBank din demonstrație și identifică adnotarea CDS. Aceasta este prima aplicare concretă, după ce studenții știu ce caută.*

## 10 Traducerea și codul genetic

O proteină este alcătuită dintr-un lanț sau mai multe lanțuri de aminoacizi. Secvența aminoacizilor contribuie la pliere și la proprietățile structurale și funcționale. Proteinele pot avea roluri enzimatice, structurale, de transport, semnalizare sau reglare.

Ribozomul citește ARN-ul mesager în grupuri de câte trei nucleotide. Fiecare grup se numește codon. Cu patru simboluri și trei poziții, avem 4 la puterea 3, adică 64 de combinații. În codul genetic standard, 61 specifică aminoacizi și trei sunt codoni stop.

Mai mulți codoni pot specifica același aminoacid. Spunem că acest cod este degenerat sau redundant; nu înseamnă că același codon are arbitrar mai multe semnificații în condițiile aceluiași cod.

Un exemplu artificial complet pentru ideea de traducere:

```text
ADN codant:  5′ ATG GAA TTT TAA 3′
ADN matriță: 3′ TAC CTT AAA ATT 5′
ARN mesager: 5′ AUG GAA UUU UAA 3′
Traducere:      Met Glu Phe Stop
Cod de o literă: M E F *
```

Stop nu este un aminoacid adăugat la proteină. Simbolul `*` îl marchează în unele rezultate informatice. AUG codifică metionina și poate funcționa ca start, dar identificarea unui început biologic de traducere depinde de context; nu orice ATG dintr-un șir este începutul unei gene.

Cadrul de citire precizează cum grupăm bazele. Dacă pornim cu o poziție mai târziu, tripletele se schimbă. Pentru ADN bicatenar există trei cadre pe fiecare orientare, ceea ce explică cele șase cadre discutate uneori în analiza secvențelor.

Nu vom traduce toate aceste cadre astăzi. Vreau să puteți explica de ce lungimea unui șir și locul de început contează, de ce apar stop-uri și de ce un transcript întreg nu este echivalent cu o CDS validată.

*Exercițiu: GAA și GAG codifică ambele glutamat în codul standard. O schimbare A → G la a treia poziție modifică ADN-ul, dar modifică neapărat aminoacidul?*

## 11 Variație genetică și mutații

Secvențele diferă între indivizi și între specii. Pot apărea substituții, în care o bază este înlocuită, inserții, în care apare material suplimentar, și deleții, în care lipsește un segment față de secvența comparată. Există și modificări structurale mai mari, precum duplicări și inversii.

Termenul variantă descrie o diferență de secvență, de obicei în raport cu o referință sau cu alte secvențe. Mutație poate desemna schimbarea sau procesul care a produs-o. Niciunul dintre termeni nu garantează, prin el însuși, un efect negativ. O variantă poate fi neutră, poate influența o funcție sau poate avea un efect dependent de context.

Într-o regiune codantă, o substituție poate fi sinonimă dacă păstrează aminoacidul, missense dacă îl schimbă, sau nonsense dacă produce un stop prematur. Un efect sinonim asupra aminoacidului nu demonstrează absența oricărui efect biologic, deoarece pot exista alte consecințe, de exemplu asupra procesării ARN-ului.

O inserție sau deleție în CDS poate deplasa cadrul de citire dacă lungimea ei nu este multiplu de trei. Acest efect se numește frameshift. Dacă lungimea este multiplu de trei, cadrul se păstrează după modificare, dar asta nu garantează că proteina funcționează normal.

SNP este un termen folosit pentru variații la un singur nucleotid în populație; SNV descrie mai general o variantă la un singur nucleotid. În practică terminologia poate varia între surse, iar baza dbSNP conține și alte tipuri de variante mici. Nu deducem tipul exact al unui record doar din numele bazei de date.

Diferențele observate într-un fișier pot proveni și din erori de secvențiere sau procesare. Alinierea ne ajută să le localizăm; validarea unei variante necesită dovezi suplimentare.

*Întrebare: dacă într-o comparație apare un gap, ce putem spune sigur și ce nu putem spune încă despre originea biologică a diferenței?*

## 12 Genetică genomică și celelalte niveluri de analiză

Genetica studiază genele, variația și transmiterea caracteristicilor. Genomica privește genomul la scară largă, inclusiv organizarea, funcția și relațiile dintre elementele lui. Domeniile se suprapun; diferența este adesea de scară și de întrebări. [NHGRI despre genomică](https://www.genome.gov/about-genomics/fact-sheets/A-Brief-Guide-to-Genomics)

Genomica structurală urmărește organizarea genomului. Genomica funcțională investighează rolurile elementelor și interacțiunile lor. Genomica comparativă examinează asemănări și diferențe între genomuri. Genomica translațională urmărește folosirea cunoștințelor în aplicații medicale. Termenul „translațională” din acest context nu trebuie confundat cu traducerea ARN-ului în proteină.

Transcriptomul reprezintă ansamblul transcripturilor dintr-un context biologic. Proteomul descrie ansamblul proteinelor dintr-un context. Aceste niveluri depind de tipul celular, condiție și moment. De aceea, o matrice de expresie trebuie însoțită de informații despre probe și măsurare.

Cantitatea de ARN nu determină direct și perfect cantitatea proteinei: intervin traducerea, degradarea și alte mecanisme. O măsurare la un nivel oferă o perspectivă, nu o descriere completă a celulei.

Bioinformatica folosește metode informatice și statistice pentru organizarea și analiza acestor date. Putem lucra cu secvențe, tabele de expresie, rețele de interacțiuni, structuri sau adnotări. Fiecare tip de date are operații și presupuneri proprii.

Astăzi rămânem în principal la secvențe. Această bază ne va ajuta ulterior să înțelegem de ce într-un laborator avem FASTA, într-altul FASTQ și într-altul o matrice cu gene pe linii și probe pe coloane.

*Revino la întrebarea despre neuron și celula musculară. Studenții ar trebui acum să poată propune comparația transcriptomului sau proteomului și să precizeze ce informație ar aduce.*

## 13 Cum ajunge biologia într-un fișier

Secvențierea este procesul experimental prin care determinăm ordinea bazelor. Multe tehnologii produc citiri, numite reads, din fragmente ale moleculelor prezente în probă. O citire nu este automat o genă și nu reprezintă neapărat un cromozom complet.

Prin asamblare putem reconstrui secvențe mai lungi din citiri, iar prin mapare putem identifica unde se potrivesc citirile pe o referință. Adnotarea adaugă interpretări: unde sunt genele, CDS-urile sau alte elemente. Acești pași nu sunt sinonimi.

O secvență de referință este un reper, nu un „genom perfect” care epuizează variația unei specii. Accesia identifică un record, iar versiunea identifică o anumită versiune a secvenței. Pentru reproducere, `NM_000546.6` este mai precis decât un simplu nume de genă.

FASTA păstrează un antet și secvența. FASTQ păstrează și scoruri de calitate pentru baze, ceea ce devine relevant la citiri. GenBank poate include secvența și adnotările ei. PDB este asociat cu structuri tridimensionale; nu este echivalent cu un simplu FASTA proteic.

```text
>exemplu_artificial descriere scurtă
ATGGAATTTTAA
```

Într-un multi-FASTA avem mai multe înregistrări. O secvență poate fi împărțită pe mai multe linii pentru lizibilitate. Începutul unei noi înregistrări este marcat de antet, nu de fiecare linie nouă.

Înainte de analiză verificăm organismul, molecula, accesie și versiune, lungimea, simbolurile și proveniența. Nu folosim numele fișierului drept singura dovadă. O secvență de transcript poate fi reprezentată în baza de date cu T, deși molecula descrisă este ARN.

*Abia aici rulează demonstrația BRCA1: identificator NM_007294.4, transcript, 7088 nucleotide, adnotare CDS. Numele genei servește drept exemplu de record, nu ca temă centrală a lecției.*

## 14 De la secvență la aliniere

Acum putem formula problema informatică în termeni biologici. Avem două secvențe și vrem să propunem poziții corespunzătoare între ele, păstrând ordinea caracterelor. Pot exista substituții, inserții sau deleții, deci nu comparăm automat indicele i din primul fișier cu indicele i din al doilea.

```text
Secvența A: A C G T
Secvența B: A - G T
```

Gap-ul ajută la reprezentarea corespondenței. Nu este o bază măsurată. Alinierea este o ipoteză despre corespondențe, evaluată printr-un sistem de scor. Un algoritm poate găsi soluția optimă pentru un model dat, însă asta nu demonstrează singur că fiecare coloană reflectă exact istoria evolutivă.

Două secvențe omoloage au o origine evolutivă comună. Ortologii sunt legați de evenimente de separare a speciilor, iar paralogii de duplicarea genelor. Similaritatea observată poate susține o ipoteză de omologie, dar homologia nu este un procent: procentul poate descrie identitatea de secvență.

La proteine, unele regiuni formează domenii structurale sau funcționale. Un motiv este un tipar de secvență, de regulă mai scurt; termenii nu sunt interschimbabili. Regiunile conservate pot orienta investigația funcției, fără să demonstreze singure rolul exact.

De aici putem continua cu slide-urile despre aliniere globală, locală, scoruri și programare dinamică. Întâi folosim secvențele artificiale `TTTACGTAAA` și `GGGACGTCCC`, ca să vedem clar mecanismul. Apoi putem compara fragmente biologice sau un set de transcripturi TP53, cu identitățile verificate.

BRCA1 și TP53 sunt exemple distincte în materialele practice. Nu le prezentăm ca pe o pereche de ortologi și nu promitem o aliniere biologic relevantă între ele doar pentru că ambele sunt gene umane.

## 15 Verificare înainte de partea practică

Înainte de a rula algoritmi, cere studenților să explice cu propriile cuvinte următoarele. Răspunsurile de mai jos sunt repere pentru discuție, nu formulări de memorat.

| Întrebare | Ideea care trebuie să apară |
| --- | --- |
| Care este relația dintre nucleotid, genă și cromozom? | Nucleotidele alcătuiesc ADN-ul; gena este o regiune; cromozomul organizează o moleculă lungă de ADN cu proteine |
| Genomul conține doar gene codante? | Include și numeroase regiuni necodante |
| De ce celule cu genom apropiat au funcții diferite? | Expresie și reglare diferite, context celular |
| Replicarea și transcripția fac același lucru? | Prima copiază ADN; a doua produce ARN dintr-o matriță ADN |
| Orice ARN este tradus? | Există ARN-uri funcționale necodante |
| Orice exon este integral CDS? | Exonii pot include regiuni netraduse |
| Alela și izoforma sunt același lucru? | Versiune de secvență la locus versus forme de transcript/proteină |
| O substituție modifică obligatoriu aminoacidul? | Codul genetic are redundanță; contează și regiunea și cadrul |
| Aceeași fracție GC implică aceeași secvență? | Compoziția nu determină ordinea |
| Două secvențe de aceeași lungime sunt deja aliniate? | Lungimea egală nu stabilește corespondențele biologice |
| O variantă față de referință înseamnă boală? | Efectul necesită dovezi și context |
| Ce păstrăm pentru reproducere? | Accesii și versiuni, surse, fișiere, parametri și comenzi |

## Cum se integrează cu prezentarea existentă

La slide-ul 1 anunță traseul molecule → informație → date → analiză. La slide-ul 2 explică obiectivele și diagnostichează vocabularul clasei. În jurul slide-ului 3 parcurge celula, ADN-ul, orientarea, cromozomii, genele, genomul și distincția genetică/genomică. Slide-ul 4 merită timp pentru replicare versus expresie, transcripție, splicing, CDS, traducere și variație. Slide-ul 5 leagă aceste noțiuni de fișiere și de primele comenzi. Slide-urile 6–15 dezvoltă alinierea pe această bază.

Aceste detalii depășesc conținutul scris pe primele slide-uri. Folosește tabla sau 4–6 slide-uri suplimentare cu structura ADN, nivelurile de organizare, splicing și exemplul de traducere. Nu înghesui tot textul în prezentare. La fiecare bloc alternează explicația cu o întrebare sau cu un exemplu mic.

Pentru o sesiune de aproximativ 110 minute, o variantă realistă este: 40–45 minute fundamente, 10 minute fișiere și verificarea mediului, 20 minute tipuri de aliniere și scoruri, 20 minute matrice și demonstrație ghidată, 10–15 minute exercițiu și recapitulare. Dacă studenții pornesc de la zero, implementarea completă a ambilor algoritmi, MSA și accesul live la baze de date vor necesita continuare. Prioritizează înțelegerea relațiilor și un exemplu corect de aliniere în această întâlnire.

## Lecturi pentru verificare și aprofundare

- [NHGRI Talking Glossary](https://www.genome.gov/genetics-glossary) pentru vocabular și ilustrații.
- [NCBI Bookshelf Molecular Biology of the Cell](https://www.ncbi.nlm.nih.gov/books/NBK21054/) pentru mecanismele moleculare.
- [EMBL-EBI despre alinierea secvențelor](https://www.ebi.ac.uk/training/online/courses/guide-to-sequence-analysis-tools/sequence-alignment/pairwise-sequence-alignment/) pentru trecerea la algoritmi.
- [Biopython Tutorial](https://biopython.org/docs/latest/Tutorial/index.html) pentru reprezentări și operații informatice.

Acest primer urmărește înțelegerea datelor și nu oferă interpretări clinice ale variantelor.
