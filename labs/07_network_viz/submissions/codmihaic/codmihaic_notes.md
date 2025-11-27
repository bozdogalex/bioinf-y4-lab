Demo01:
Hub genes: [('GeneC', 3), ('GeneA', 2)]
Am salvat demo_network.png și demo_hubs.csv

Ex01:
Grafic: 998 noduri, 3070 muchii.
Am salvat figura în: labs/07_networkviz/submissions/codmihaic/network_codmihaic.png
Am salvat hub genes în: labs/07_networkviz/submissions/codmihaic/hubs_codmihaic.csv

1) Metoda de layout folosită:
În script am folosit nx.spring_layout, adică un layout de tip spring / force-directed (Fruchterman–Reingold).
Pe scurt, algoritmul tratează muchiile ca niște „arcuri” care trag nodurile conectate unul spre celălalt și introduce forțe de respingere între toate nodurile. Rezultatul este o poziționare în plan 2D în care:
- genele puternic conectate ajung apropiate spațial (formează „insule”);
- genele cu puține conexiuni sunt împinse spre marginea grafului;
- layout-ul este relativ intuitiv vizual pentru rețele biologice de mărime medie.
Am folosit și un număr fix ca seed (42) pentru a obține același layout de fiecare dată, ceea ce ajută la compararea și reproducerea rezultatelor.

2) Reflecție: ce avantaje aduce vizualizarea față de analiza numerică din Lab 6?
În Lab 6 am lucrat mai ales numeric: matrici de corelație, praguri, matrice de adiacență, liste de module (gene → modul), eventual scoruri de centralitate. Aceste rezultate sunt corecte, dar greu de „intuit” doar din tabele, pe când vizualizarea rețelei aduce câteva avantaje clare:
- înțelegere globală mai rapidă: dintr-o singură figură văd imediat câte module mari/mici există, cât de bine sunt separate și dacă există module „aglomerate” sau foarte dispersate.
- localizarea hub-urilor și a legăturilor-cheie: genele hub apar central, cu multe muchii; pot observa și gene punte între module (posibile integratoare de procese biologice), lucru greu de dedus doar din liste de grade sau corelații.
- detectarea artefactelor sau a pattern-urilor neașteptate: dacă un modul apare vizual „lipit” de altul sau dacă există componente izolate, pot suspecta probleme de prag, de filtrare sau pot formula ipoteze biologice noi.
- comunicare mai ușoară: o figură cu module colorate și hub-uri etichetate este mult mai ușor de prezentat și explicat decât o matrice numerică sau un fișier CSV.
Pe scurt, analiza numerică din Lab 6 ne spune „cât de puternice sunt conexiunile și cum sunt grupate genele”, iar vizualizarea din Lab 7 ne arată „cum arată efectiv rețeaua”, făcând mult mai ușoară interpretarea biologică și verificarea intuitivă a rezultatelor.