### Raport Vizualizare Retea (Lab 7)

**Metoda folosita**

Pentru vizualizare am folosit functia spring_layout din NetworkX. Asta ne aranjeaza nodurile (genele) in spatiu in asa fel incat genele care se "trag" (adica sunt conectate prin corelatie mare) sa stea aproape. Nodurile au fost colorate strict pe baza modulului identificat in **Lab 6 (WGCNA)**. Am marit vizibil nodurile care au un grad mare (Hub Genes) ca sa le vad repede pe cele mai conectate din retea.

---

**Avantaje vizualizare vs. analiza numerica**

Vizualizarea retelei face interpretarea mult mai usoara decat simpla analiza numerica a matricilor. In Lab 6 am avut tabele, dar nu vedeam structura. Acum, vizualizarea arata imediat **densitatea** modulelor si **conexiunile cheie**.

De exemplu, observ ca modulele importante au tendinta de a forma grupuri bine definite. Pot vedea dintr-o privire care module sunt mai "interconectate" decat altele. De asemenea, genele Hub (nodurile mari) sar in ochi si pot fi analizate prima data, lucru care e greu doar dintr-un CSV mare cu Grade. Vedem direct care gene sunt centrale in fiecare grup.