Precum a fost mentionat in fisierul README.md am folosit datele descarcate in laboratorul 1, acestea aflate in fisierul my_tp53.fa.


In urma rularii urmatoarelor comenzi: 
python ex01_global_nw.py --fasta ../../../my_tp53.fa --i1 0 --i2 1
python ex02_global_nw.py --fasta ../../../my_tp53.fa --i1 0 --i2 1

am tras urmatoarele concluzii: 

Aliniere globală (Needleman–Wunsch) este potrivită când vrei aliniere end-to-end (întregi secvențe) și secvențele sunt de lungimi similare și se așteaptă omologie pe aproape toată lungimea.

Aliniere locală (Smith–Waterman) este potrivită când interesezi regiuni conservate sau domenii comune între secvențe care altfel pot fi foarte diferite (ex.: un domeniu conservat în proteine din ORFs diferite, fragmente de read-uri vs genom).