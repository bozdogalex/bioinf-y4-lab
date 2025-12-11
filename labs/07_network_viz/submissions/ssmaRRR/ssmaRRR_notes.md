Am facut exercitiul, si mi-a generat urmatoarele fisiere: network_ssmaRRR.png si hubs_ssmaRRR.csv pe care le-am salvat in 
labs/07_network_viz/submissions/ssmaRRR

Am folosit metoda de layout spring, prin urmatoarea linie de cod: pos = nx.spring_layout(G, seed=42)

Vizualizarea rețelei completează analiza numerică din Lab 6, pentru că îmi permite să văd direct structura rețelei, nu doar să citesc module și grade în tabele. Gruparea spațială a nodurilor evidențiază clar modulele dense, ceea ce face mai intuitiv de înțeles organizarea globală a rețelei.