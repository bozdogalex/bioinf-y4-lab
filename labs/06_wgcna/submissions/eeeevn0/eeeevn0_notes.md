**ex01**

Pentru acest laborator am folosit un fișier de expresie derivat din setul de date **GSE115469**.Fișierul original era foarte mare, așa că am extras un subset de 3000 de gene și l-am salvat la: data/work/eeeevn0/lab06/expression_matrix.csv, ca să putem lucra cu el mult mai ușor.

**Metrica folosita**
Am folosit corelația Spearman, deoarece se potrivește bine datelor RNA-Seq și arată relațiile dintre gene chiar dacă valorile nu sunt distribuite liniar.

**Pragul ales**
Pentru construirea rețelei co-expresie am aplicat un prag de 0.6 pe valoarea absolută a corelației. Am vrut să păstrăm doar legăturile mai puternice dintre gene, așa că toate corelațiile sub 0.6 au fost eliminate din rețea.

După preprocesare (log-transformare și filtrarea genelor cu varianță foarte mică) au rămas 50 de gene. Rețeaua finală a avut 50 de noduri și 219 de muchii, iar algoritmul Louvain a detectat 6 module. De exemplu, genele ribozomale, cum ar fi RPL13, RPL19, RPS23 și RPS29, au apărut împreună în același modul, ceea ce are sens deoarece toate sunt implicate în sinteza proteinelor. La fel, genele imune precum NKG7 și CTSW au fost grupate în același modul, fiind legate de funcțiile celulelor imunitare.


**Cum diferă o rețea de co-expresie față de clustering-ul clasic?**  
Clustering-ul grupează genele doar după cât de mult se aseamănă între ele, privindu-le ca pe un singur întreg, în timp ce o rețea de co-expresie arată legăturile directe dintre gene, adică ce gene se corelează puternic între ele iar modulele se formează din aceste conexiuni.
