# Note de predare pentru laboratoarele 1 și 2 de bioinformatică

Textul urmărește cele 15 slide-uri ale prezentării „From Genomics to Sequence Alignment”. Firul lecției este înțelegerea moleculelor și a informației genetice, apoi trecerea la fișiere și algoritmi. TP53 rămâne un exemplu opțional în partea practică. Exemplele foarte scurte sunt artificiale și servesc explicării algoritmilor; secvențele din repository sunt folosite separat, pentru analiza datelor biologice.

Poți folosi paragrafele ca discurs și poți scurta explicațiile în funcție de răspunsurile studenților. Indicațiile scrise cu italice sunt pentru tine. Întrebările sunt adresate clasei. Comenzile actuale sunt în README-urile laboratoarelor 1 și 2; pentru bazele biologice folosește primerul asociat.

## Slide 1 De la biologie la date și algoritmi

Astăzi construim o legătură între două lumi. În biologie avem celule, molecule și procese. În calculator avem fișiere, șiruri de caractere și algoritmi. Ca să folosim corect instrumentele, trebuie să înțelegem ce reprezintă datele.

Să începem cu o întrebare simplă: cum se explică faptul că un neuron și o celulă musculară pot avea aproape același genom, dar funcționează atât de diferit? Răspunsul ne va conduce către diferența dintre informația genetică disponibilă și informația folosită într-un anumit context.

Nu presupun că termenii de biologie sunt deja familiari. Vom clarifica mai întâi relația dintre ADN, cromozom, genă și genom. Apoi vedem cum ajungem de la ADN la ARN și proteină și de ce există fișiere diferite pentru aceste niveluri.

Abia după ce știm ce conține un fișier discutăm compararea secvențelor. Vom începe cu exemple artificiale scurte, care pot fi urmărite pe tablă. Mai târziu folosim și date biologice, ca aplicație.

*Suportul detaliat pentru explicațiile introductive este [primerul de genetică și genomică](genetics-genomics-primer-ro.md). Folosește secțiunile lui ca discurs, nu doar lista de definiții. Exemplele BRCA1 și TP53 apar după fundamente.*

## Slide 2 Ce trebuie să putem explica la final

La final aș vrea să puteți explica trei lucruri: ce obiect biologic avem, cum este reprezentat și ce întrebare putem investiga cu el.

O secvență este ordinea unor unități: nucleotide în ADN și ARN, aminoacizi în proteine. O genă este o regiune de ADN cu informație pentru un produs funcțional. Un cromozom organizează o moleculă lungă de ADN împreună cu proteine. Genomul include întregul material genetic, cu regiuni codante și necodante.

Aceste obiecte sunt legate, dar nu sunt sinonime. Un fișier de transcript nu reprezintă automat întreaga genă genomică, iar o secvență proteică nu se compară ca și cum literele ei ar reprezenta nucleotide.

Pe parcurs trebuie să înțelegem transcripția, procesarea ARN-ului și traducerea suficient cât să evităm aceste confuzii. La partea informatică, vom recunoaște formatele, vom calcula proprietăți simple și vom interpreta o aliniere.

Nu vă cer să memorați numele tuturor enzimelor. Vă cer să puteți urmări traseul informației și să explicați ce face fiecare operație din program.

*Pe slide, „genus” înseamnă gen taxonomic; dacă intenția era „gene”, corectează termenul. Înainte să mergi mai departe, cere o propoziție despre relația genă–cromozom–genom. Dacă răspunsurile sunt neclare, folosește primerul, secțiunile 2–6.*

## Slide 3 Ce studiază genetica și genomica

Genetica studiază genele, variația și moștenirea. Genomica privește materialul genetic la scară largă: organizarea, funcția și relațiile dintre elementele lui. Genomul conține mai mult decât regiunile care codifică proteine.

Înainte să discutăm ramurile genomicii, trebuie să putem desena nivelurile de organizare. Într-o celulă eucariotă, cea mai mare parte a ADN-ului se află în nucleu. ADN-ul este organizat în cromozomi; genele sunt regiuni ale acestui ADN. Există și ADN în organite, cum sunt mitocondriile.

ADN-ul este alcătuit din nucleotide. Literele A, C, G și T reprezintă bazele. Catenele sunt complementare și antiparalele, iar secvențele sunt de regulă raportate în direcția 5′ spre 3′. Orientarea contează în comparațiile informatice.

Acum putem înțelege ramurile de pe slide. Genomica structurală privește organizarea, cea funcțională urmărește rolurile elementelor, cea comparativă examinează asemănările și diferențele, iar cea translațională urmărește aplicațiile medicale.

Datele provin din măsurători experimentale și din prelucrarea lor. Înainte de a aplica un algoritm, verificăm ce tip de secvență avem și din ce context provine. Această verificare este o parte a analizei.

*Dezvoltă aici primerul 2–6 și 12. Arată pe tablă o regiune de cromozom cu mai multe gene și regiuni dintre ele, apoi o dublă catenă scurtă. Explică alela ca versiune la un locus; nu o confunda cu o genă diferită.*

## Slide 4 Cum este utilizată informația genetică

Schema ADN → ARN → proteină descrie traseul clasic pentru genele codante. Pentru a o înțelege, trebuie să separăm replicarea de expresie.

Replicarea copiază ADN-ul, de exemplu înaintea diviziunii. Transcripția produce ARN folosind o catenă de ADN ca matriță. Traducerea folosește informația din ARN-ul mesager pentru a construi un lanț de aminoacizi.

ARN-ul nu este doar o etapă intermediară. Există ARN-uri ribozomale, de transfer și de reglare, cu propriile funcții. Nici expresia genică nu înseamnă că întregul genom este folosit uniform: activitatea diferă între celule și condiții.

La multe gene eucariote, ARN-ul este procesat. Intronii sunt eliminați, exonii sunt uniți, iar transcriptul matur poate conține UTR-uri și o regiune codantă CDS. Un exon nu este obligatoriu integral codant. Splicingul alternativ poate produce transcripturi diferite.

La traducere citim câte trei nucleotide, adică un codon. Exemplul artificial ATG GAA TTT TAA, dacă reprezintă ADN codant în cadrul potrivit, corespunde ARN-ului AUG GAA UUU UAA și produsului Met–Glu–Phe urmat de stop.

O substituție poate păstra aminoacidul, îl poate schimba sau poate produce stop. Un indel poate deplasa cadrul de citire. De aceea, o diferență de secvență nu are aceeași semnificație în orice regiune.

*Acest slide merită cel mai mult timp teoretic. Folosește primerul 7–11, cu desenul splicingului și exemplul de traducere. Rulează apoi demo02_seq_ops.py și explică stop-urile interne din exemplul lui. Programul operează pe șirul oferit; nu identifică automat gene, CDS-uri sau cadre biologic corecte.*

## Slide 5 De la molecule la fișiere

După aceste explicații, putem deschide primul fișier cu o întrebare precisă: ce obiect reprezintă?

FASTA conține un antet care începe cu semnul mai mare și o secvență. Într-un multi-FASTA avem mai multe înregistrări. Secvența poate fi împărțită pe mai multe linii; antetul marchează începutul unei alte înregistrări.

FASTQ include și calitatea identificării bazelor. GenBank poate conține secvența împreună cu adnotări, inclusiv CDS. PDB este asociat structurii tridimensionale. Alegerea formatului depinde de ce informație vrem să păstrăm.

Înainte de calcul verificăm organismul, tipul moleculei, accesie și versiune, lungimea și simbolurile. Acum putem folosi BRCA1 drept exemplu de transcript real: recordul inclus în demonstrație are accesie NM_007294.4 și 7088 nucleotide. Numele genei este un exemplu de metadată; principiile se aplică și altor gene.

Putem calcula fracția GC, dar aceasta descrie compoziția. Două secvențe cu aceeași compoziție pot avea ordine diferită. Dacă vrem să comparăm poziții și diferențe, avem nevoie de aliniere.

*Rulează demo-ul BRCA1 fără --refresh pentru lectura offline. Cere identificarea moleculei și a CDS-ului în GenBank. Abia acum introdu opțional fișierele TP53 ca alte exemple de secvențe. Nu prezenta BRCA1 și TP53 drept ortologi. Continuă cu exemple artificiale pentru algoritmi.*


## Slide 6 De ce avem nevoie de aliniere

Să presupunem că am verificat datele și avem secvențele corecte. Vrem să testăm ipoteza despre asemănarea lor. De ce nu comparăm pur și simplu primul caracter cu primul, al doilea cu al doilea și așa mai departe?

Pentru că secvențele se pot modifica și prin inserții și deleții. Dacă într-o secvență apare un caracter suplimentar, toate caracterele de după el se pot deplasa cu o poziție. O comparație rigidă ar putea transforma un singur eveniment într-un șir lung de nepotriviri.

Alinierea propune o corespondență între poziții. Această corespondență ne permite să discutăm diferențele într-un mod mai util biologic.

Pe slide avem câteva aplicații. Putem căuta indicii de origine evolutivă comună, adică omologie. Putem observa substituții și regiuni compatibile cu inserții sau deleții. Putem formula ipoteze despre funcția unei secvențe necunoscute, prin comparație cu una caracterizată. Putem pregăti datele pentru analize filogenetice.

În exemplul TP53, ne interesează atât asemănarea generală, cât și regiunile care se păstrează între organisme. O regiune conservată poate fi un punct bun de pornire pentru o întrebare despre funcție.

Dar asemănarea observată trebuie interpretată în context: lungime, compoziție, acoperire și metoda folosită. O potrivire foarte scurtă poate apărea întâmplător.

Și o precizare de vocabular: putem măsura un procent de identitate sau de similaritate. Omologia descrie originea comună; nu spunem că două secvențe sunt „80% omoloage”.

*Pe tablă: compară ACGT cu AGT, mai întâi fără gap și apoi cu A-GT. Cere-le să explice ce s-a schimbat în interpretare.*

## Slide 7 Ce face efectiv alinierea

În exemplul de pe tablă, am obținut:

```text
A C G T
A - G T
```

Am păstrat ordinea caracterelor din fiecare secvență și am introdus un gap pentru a propune poziții corespunzătoare. Gap-ul este un simbol al alinierii. Nu este o bază pe care am măsurat-o în probă.

Comparația este compatibilă cu o inserție sau o deleție. Ca să stabilim direcția evolutivă a schimbării, avem nevoie de informații suplimentare. Din două secvențe singure nu aflăm automat care reprezintă starea ancestrală.

Pentru secvențe mai lungi, pot exista foarte multe moduri de a introduce gap-uri. Unele alinieri arată plauzibil, altele au atât de multe întreruperi încât devin greu de interpretat. Avem nevoie de o regulă comună pentru a le compara.

Această regulă este funcția de scor. Definim cât valorează o potrivire, cât costă o nepotrivire și cât costă un gap. Algoritmul caută apoi o aliniere care maximizează scorul.

Observați că i-am dat algoritmului un model al problemei. Dacă schimbăm regulile, se poate schimba și soluția. Prin „aliniere optimă” înțelegem optimă pentru secvențele, tipul de aliniere și scorurile alese. Pot exista mai multe soluții cu același scor maxim.

De aceea, când comparați rezultatul vostru cu al unui coleg, începeți prin a verifica parametrii. Două rezultate diferite nu înseamnă neapărat că unul dintre programe este greșit.

*Întrebare: de ce nu recompensăm potrivirile și lăsăm toate gap-urile gratuite? Lasă-i să observe că am putea favoriza alinieri fragmentate, cu multe întreruperi.*

## Slide 8 Alegerea între global local și semi global

Înainte să rulăm programul, trebuie să precizăm ce fel de răspuns căutăm.

Dacă avem două secvențe comparabile ca lungime și ne așteptăm să fie înrudite de-a lungul întregii secvențe, alegem o aliniere globală. Aceasta include ambele secvențe integral. Needleman–Wunsch este algoritmul clasic pentru această problemă.

Dacă presupunem că doar o regiune este comună, alinierea locală este mai potrivită. Putem avea, de exemplu, două proteine care împart un domeniu, dar diferă în rest. Smith–Waterman caută perechea de subsecvențe cu scorul maxim.

Alinierea semi-globală este utilă în situații în care anumite porțiuni de la capete nu trebuie penalizate. Dacă potrivim un fragment într-o referință mai lungă, vrem să putem ignora flancurile referinței care nu fac parte din fragment. Există mai multe configurații ale capetelor libere; spunem explicit ce lăsăm nepenalizat.

Să aplicăm ideea la datele noastre. Două transcripturi pot conține regiuni netraduse de lungimi diferite. O aliniere globală a transcripturilor complete și o aliniere a regiunilor codante răspund unor întrebări diferite. Alegerea datelor și alegerea metodei trebuie făcute împreună.

Nu trebuie să memorăm „global este bun” sau „local este bun”. Trebuie să putem spune de ce folosim o anumită metodă pentru problema dată.

*Cere trei răspunsuri rapide: două versiuni complete ale aceleiași gene; un domeniu comun între proteine altfel diferite; un fragment scurt într-o referință lungă. Discută și presupunerile necesare, nu doar etichetele.*

## Slide 9 Cum definim scorul

Am stabilit ce vrem să aliniem. Acum trebuie să stabilim cum evaluăm o aliniere.

Pentru un exemplu de ADN, putem folosi plus unu pentru o potrivire, minus unu pentru o nepotrivire și minus doi pentru fiecare poziție de gap. Alegem aceste valori pentru a putea urmări calculul ușor.

```text
A C G T
A - G T
```

Avem trei potriviri și un gap, deci scorul este trei minus doi, adică unu. Dacă modificăm costul gap-ului, se modifică scorul acestei alinieri. În alte exemple se poate modifica și alinierea care câștigă.

La proteine, folosirea aceleiași penalizări pentru toate perechile de aminoacizi diferiți ar pierde informație. Unele înlocuiri sunt mai compatibile cu proprietățile și rolul regiunii respective decât altele. Matricile BLOSUM și PAM oferă scoruri pentru perechi de aminoacizi, pe baza modelelor și observațiilor din care au fost construite.

Pentru gap-uri distingem penalizarea liniară de cea afină. În modelul liniar, fiecare poziție are același cost. În modelul afin, avem un cost pentru deschiderea gap-ului și unul pentru extindere. O convenție uzuală pentru un gap de lungime k este costul de deschidere plus k minus unu ori costul de extindere; trebuie verificată convenția instrumentului folosit.

Intuiția este că un singur fragment inserat poate fi tratat diferit de mai multe inserții separate. Modelul afin permite această diferență.

În exercițiile noastre vom începe cu modelul liniar. Când comparați implementarea proprie cu Biopython, folosiți aceleași scoruri. În repository, valorile implicite ale demonstrației, ale exercițiului global și ale exercițiului local diferă, deci scorurile afișate nu trebuie comparate direct.

*Corecție de slide: BLOSUM și PAM sunt denumirile corecte. Referință pentru scoruri și tipuri de aliniere: [EMBL-EBI](https://www.ebi.ac.uk/training/online/courses/guide-to-sequence-analysis-tools/sequence-alignment/pairwise-sequence-alignment/).*

## Slide 10 Cum găsește algoritmul alinierea

Pentru două secvențe de câteva litere putem propune alinieri manual. Pentru secvențe lungi, numărul posibilităților crește foarte mult. Programarea dinamică ne permite să construim soluția folosind rezultate pentru probleme mai mici.

Punem o secvență pe linii și cealaltă pe coloane. În fiecare celulă păstrăm cel mai bun scor pentru prefixele corespunzătoare. Asta înseamnă că o celulă nu descrie doar cele două caractere aflate în dreptul ei; rezumă cea mai bună comparație până în acel punct.

Pentru a o calcula avem trei posibilități. Din diagonală, aliniem câte un caracter din fiecare secvență. De sus, consumăm un caracter din secvența de pe linii și îl aliniem cu un gap. Din stânga, consumăm un caracter din secvența de pe coloane și îl aliniem cu un gap. Alegem varianta cu scorul maxim.

Haideți să facem o matrice foarte mică pentru AC și AG, cu match plus unu, mismatch minus unu și gap minus doi. Prima linie și prima coloană sunt zero, minus doi, minus patru. Celula A cu A primește unu. Pentru C cu G, varianta diagonală dă unu minus unu, adică zero; în acest exemplu este cea mai bună variantă. Scorul final este zero.

După completare, facem traceback: refacem traseul care a produs scorul. Pentru global pornim din colțul din dreapta jos. Pentru local inițializăm cu zero, permitem și zero în recurență, apoi pornim din celula maximă și ne oprim la zero.

Cu lungimi n și m, avem aproximativ n ori m celule de calculat. La șapte litere este foarte puțin; la mii de nucleotide diferența se simte, mai ales în implementarea didactică Python.

*Pe tablă, matricea globală completă este: prima linie 0, −2, −4; a doua −2, 1, −1; a treia −4, −1, 0. Arată o singură celulă în detaliu și cere clasei următoarea. Recursia simplă de pe slide folosește gap-uri liniare.*

## Slide 11 Ce vedem diferit în rezultatele globale și locale

Vreau acum să anticipăm rezultatul înainte să îl calculăm. Folosim două secvențe artificiale, construite ca să aibă o regiune comună și flancuri diferite:

```text
TTTACGTAAA
GGGACGTCCC
```

Ce se va întâmpla la alinierea globală? Ambele secvențe trebuie incluse integral. Pot apărea nepotriviri și, în funcție de scoruri, gap-uri. Ce se va întâmpla la alinierea locală? Algoritmul poate alege regiunea ACGT și lăsa în afara alinierii capetele care reduc scorul.

Un rezultat local scurt, dar foarte bun, nu ne spune automat că secvențele complete sunt foarte asemănătoare. Ne spune că există o regiune cu un scor bun. De aceea, când prezentați rezultatul, menționați și cât din secvență a fost aliniat.

Acum mărim problema. Dacă avem o secvență și vrem să căutăm asemănări într-o bază de date foarte mare, comparația exactă cu fiecare intrare poate fi costisitoare. BLAST folosește o strategie euristică: caută potriviri scurte și extinde regiunile promițătoare. Asta îl face practic pentru căutări la scară mare. [Documentația NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-topics/blastsearchparams.html)

Smith–Waterman și BLAST nu sunt același algoritm. Primul rezolvă problema exactă a alinierii locale pentru modelul dat; BLAST folosește euristici pentru căutare rapidă și poate rata unele potriviri.

*Rulează demo01_pairwise.py pe perechea artificială pregătită în ghid. Apoi rulează pe primele 60 de nucleotide din fișierul TP53. Precizează că acel prefix este un exemplu tehnic; el nu reprezintă automat o regiune codantă sau un domeniu funcțional.*

## Slide 12 Ce aduce a treia secvență

Am comparat două secvențe și am identificat o regiune comună. Ce câștigăm dacă adăugăm o a treia secvență?

O aliniere multiplă ne permite să privim aceeași coloană în mai multe organisme. Unele poziții pot varia frecvent, altele pot rămâne identice. Regiunile conservate ne pot orienta spre elemente importante pentru funcție sau structură.

În exemplul cu p53, ne interesează dacă există segmente păstrate între organisme. Totuși, simpla existență a unui segment identic nu demonstrează singură ce funcție are. Pentru această interpretare consultăm și adnotările, structurile sau alte rezultate experimentale.

În mod obișnuit, numim aliniere multiplă comparația a trei sau mai multe secvențe. Problema este mai dificilă decât cea pe perechi, așa că instrumentele practice folosesc frecvent euristici. Clustal Omega construiește alinierea progresiv, ghidată de relațiile estimate dintre secvențe. T-Coffee urmărește și consistența corespondențelor dintre comparații.

Arborele ghid folosit la construirea unei alinieri are un rol algoritmic. Nu îl tratăm automat ca pe o filogenie validată.

În rezultat vă cer să alegeți o regiune și să explicați ce observați: câte poziții sunt identice, unde apar gap-uri și dacă asemănarea se întinde pe mai multe coloane. Abia după observație formulăm interpretarea.

*Folosește fișierul cu cele trei transcripturi pentru aliniere de nucleotide. Setul proteic corectat conține acum p53 umană, de șoarece și de pește zebră; consultă proveniența. Salvează alinierea descărcată, fără să depinzi de un link temporar.*

Întrebare pentru voi: dacă o regiune este conservată în toate cele trei secvențe, ce ați verifica înainte să spuneți că este esențială pentru funcție?

## Slide 13 Cum organizăm analiza cu instrumentele disponibile

Să punem instrumentele în ordinea în care ne ajută la lucru.

Mai întâi avem nevoie de date și de identificatorii lor. Le putem consulta în baze de date precum NCBI și le putem citi în Python cu Biopython. Apoi alegem metoda de comparație.

Pentru alinierea a două secvențe există instrumente precum needle pentru global și water pentru local. Pentru căutări de nucleotide în baze de date avem blastn. Pentru aliniere multiplă avem Clustal Omega. Acestea sunt programe distincte; instalarea Biopython nu le instalează automat.

În repository, demonstrația de astăzi este un script Python. Îl rulăm din terminalul Codespaces, pornind din rădăcina proiectului. Codespaces ne oferă un calculator Linux în cloud, cu un mediu configurat pentru curs. Browserul este interfața prin care lucrăm.

Deschidem scriptul, identificăm fișierul de intrare, apoi rulăm comanda. După rezultat, schimbăm un parametru sau setul de date și explicăm diferența. Așa transformăm o execuție într-un experiment.

Codul actual folosește pairwise2. Acesta este depreciat, deși funcționează în mediul verificat pentru curs. Pentru implementări noi, Biopython recomandă PairwiseAligner. Avertismentul de depreciere semnalează o schimbare recomandată a API-ului; nu înseamnă că rezultatul rularii a eșuat. [Documentația Biopython](https://biopython.org/DIST/docs/tutorial/Tutorial-1.80.html)

AlignIO ne ajută să citim și să scriem alinieri. Nu trebuie confundat cu un algoritm care calculează automat orice aliniere.

La final, păstrăm codul, parametrii, identificatorii secvențelor și observațiile. Un coleg trebuie să poată reproduce rezultatul fără să ghicească ce am făcut.

*În clasa combinată, terminalul este suficient. Jupyter poate rămâne opțional. Nu introduce MLflow în această sesiune doar pentru că este instalat în imagine.*

## Slide 14 Cum evaluăm rezultatul

Avem acum o aliniere. Ce anume ne convinge că rezultatul este relevant?

Primul indicator este scorul. El rezumă alinierea în sistemul de evaluare ales. Un scor de 50 nu are o interpretare universală. Trebuie să știm scorurile pentru potriviri și nepotriviri, penalizările de gap și, la proteine, matricea de substituție.

Al doilea indicator este procentul de identitate. Dacă avem 80 de poziții identice într-o aliniere de 100 de coloane, rezultatul este 80%, folosind aceste coloane ca numitor. Dar trebuie să precizăm convenția pentru gap-uri și poziții ambigue.

Comparați două rezultate ipotetice: 100% identitate pe opt poziții și 85% identitate pe 300 de poziții. Care oferă mai multă informație despre relația dintre secvențe? Procentul singur nu este suficient. Privim lungimea și acoperirea, iar pentru o evaluare statistică folosim indicatorii potriviți.

În BLAST, E-value este numărul de rezultate cu scor cel puțin la fel de bun pe care ne așteptăm să le găsim întâmplător într-o căutare comparabilă. Un E-value mic indică o semnificație statistică mai mare. Dimensiunea bazei de date influențează acest indicator. El nu reprezintă probabilitatea ca două proteine să aibă aceeași funcție. [Explicația NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&DOC_TYPE=FAQ&PAGE_TYPE=BlastDocs)

Pentru alinierea multiplă putem calcula și un consens, alegând, într-o variantă simplă, caracterul cel mai frecvent din fiecare coloană. Consensul rezumă setul de date și nu trebuie să coincidă cu o secvență existentă.

Mai avem distanțele: dacă numărăm diferențe poziție cu poziție, trebuie să comparăm poziții corespunzătoare. Trunchierea a două secvențe la aceeași lungime nu realizează alinierea lor. Acesta este un punct pe care îl vom verifica în demonstrație.

*Întrebare: două fișiere au aceeași fracție GC, dar alinierea lor este slabă. Există o contradicție? Răspunsul ar trebui să distingă compoziția de ordinea caracterelor.*

## Slide 15 Ce putem concluziona și ce urmează

Să revenim la traseul cu care am început: cum trecem de la molecule și procese biologice la date și la o analiză informatică?

Acum avem un traseu de lucru. Identificăm organismul și tipul moleculei, verificăm secvențele și lungimile, alegem o metodă, precizăm scorurile și interpretăm alinierea. Putem spune ce regiuni se potrivesc și unde apar diferențe. Putem formula ipoteze despre conservare și despre relațiile dintre secvențe.

Aceste idei apar în multe aplicații. Compararea cu o referință ajută la identificarea variantelor candidate. Alinierile multiple oferă poziții comparabile pentru analiza filogenetică. Asemănările dintre proteine pot contribui la identificarea domeniilor și la investigarea unor funcții comune.

În metagenomică, comparăm secvențe provenite din amestecuri de organisme cu referințe cunoscute. În cercetarea pentru repoziționarea medicamentelor, asemănările dintre ținte pot contribui la formularea unor ipoteze, împreună cu alte date. O aliniere singură nu validează efectul unui medicament.

Pentru predare, vreau să formulați o observație și o interpretare și să le țineți distincte. De exemplu: „Am obținut o regiune de lungime X cu Y poziții identice” este observația. „Această conservare sugerează că regiunea merită investigată funcțional” este interpretarea.

La următorul laborator vom lucra cu citiri NGS. Vom avea multe fragmente, calități de secvențiere și o referință. Întrebarea devine: unde se potrivește fiecare fragment și ce dovezi avem pentru o diferență față de referință?

Până atunci, verificați că rezultatul de astăzi poate fi reprodus. În notițe trebuie să existe datele folosite, comanda, parametrii și concluzia voastră, inclusiv o limitare.

*Încheiere interactivă: cere fiecărei echipe trei propoziții — „Am comparat…”, „Am observat…”, „Nu putem concluziona încă…”. Revino la distincția dintre observație și interpretare și cere exemple de noțiuni biologice care au influențat analiza.*

## Ritmul sesiunii combinate

Pentru începători, alocă aproximativ 40–45 minute fundamentelor din primer înainte de a accelera spre algoritmi. Urmează 10 minute de fișiere și mediu, 20 de minute de aliniere și scoruri, 20 de minute de programare dinamică și demonstrație, apoi 10–15 minute de exercițiu și recapitulare. O implementare completă a ambilor algoritmi și MSA pot necesita continuare.

Explică o noțiune, arată un exemplu mic, apoi cere o predicție sau o explicație. Păstrează secvențele artificiale pentru urmărirea algoritmului și secvențele biologice pentru interpretare. TP53 nu trebuie să introducă vocabular nou înainte de fundamente.

## Corecții de făcut în prezentare înainte de predare

- Clarifică „genus” versus „gene”.
- Înlocuiește „rearranging sequences” cu o formulare care spune că se păstrează ordinea caracterelor și se introduc gap-uri.
- Corectează „Blossum / PAL” în „BLOSUM / PAM”.
- Înlocuiește „BLAST = quick Smith-Waterman” cu „BLAST = căutare euristică a similarității locale”.
- Definește identitatea prin raportul dintre pozițiile identice și lungimea evaluată, cu o convenție explicită pentru gap-uri.
- Definește E-value ca număr așteptat de rezultate întâmplătoare cel puțin la fel de bune, nu ca probabilitate.
- Prezintă pairwise2 drept API folosit de materialele actuale și PairwiseAligner drept alternativa recomandată pentru cod nou.
- Spune că alinierea este o etapă de bază în multe fluxuri, precedată de alegerea și verificarea datelor; nu este obligatoriu prima operație a oricărei analize genomice.

Material de bază: [prezentarea laboratorului](https://docs.google.com/presentation/d/1cPgRk1VZVSjrxnnYIMFjvDx3ktO0-Bp2EZu7M4hRY5s/edit) și [repository-ul cursului](https://github.com/bozdogalex/bioinf-y4-lab), cu corecțiile Lab 1/2 din 5 octombrie 2026.
