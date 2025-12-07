Metoda de Layout

Pentru vizualizarea rețelei de co-expresie, am folosit algoritmul Circular Layout (nx.circular_layout).

Am ales această metodă deoarece este rapidă și eficientă din punct de vedere al resurselor (CPU/Memorie) pe grafuri mari, asigurând că vizualizarea se finalizează fără erori de tip Terminated.

Reflecție: Avantajele Vizualizării (Vizual vs. Numeric)

Vizualizarea rețelelor (față de analiza numerică din Lab 6) aduce avantaje esențiale:

1.Validarea Structurii Modulelor: Prin vizualizare (colorând nodurile după modul), putem vedea imediat dacă genele din același modul se grupează fizic sau dacă sunt împrăștiate.

2.Identificarea Genelor Hub: Vizualizarea ne arată genele hub (nodurile cu dimensiunea mărită) și modul în care acestea acționează ca punți de legătură în rețea.

3.Detectarea Tiparelor Imediate: O vizualizare transformă mii de rânduri și coloane de corelații într-o hartă funcțională, permițând identificarea rapidă a subgrupurilor și a interacțiunilor neașteptate.