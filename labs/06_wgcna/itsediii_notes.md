Am ales metoda de corelație Spearman deoarece este mai robust la distribuții non-normale și la outlieri și captează relații monotone, nu doar liniare (cum face Pearson). În contextul expresiei genice (puternic skewed + zgomot), Spearman este de obicei mai stabil.

Am folosit pragul de adiacență |cor| ≥ 0.6 pentru a crea o rețea mai curată, eliminând corelațiile slabe care pot introduce zgomot.



În clusteringul clasic (cum am făcut în Lab 5), genele sunt grupate doar pe baza asemănării generale dintre profilele lor de expresie. Practic, algoritmul caută „ce seamănă cu ce” și formează câteva grupuri mari.

O rețea de co-expresie funcționează diferit: aici ne uităm la gene ca la nodurile unui graf și la corelațiile dintre ele ca la legături. Modulele nu sunt create din distanțe globale, ci din modul în care genele sunt conectate între ele. Algoritmi precum Louvain caută direct zone dense din graf (comunități).