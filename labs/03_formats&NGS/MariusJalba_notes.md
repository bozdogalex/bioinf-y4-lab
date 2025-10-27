Am rulat demo01_fastq_qc.py:
Reads: 31
Mean length: 151.00
N rate: 0.0004
Mean Phred: 35.82

Am rulat demo02_fastq_stats.py:
Read ATGCTAGC mapped at position 0
Read GATCGATC mapped at position 19
Read TACGATCG mapped at position 28
Read GGGGGGGG did not map

Am rulat ex01_fetch_fastq.py:
am descarcat FASTQ pentru accession SRR684066
fisierul a fost descarcat in /data/work/MariusJalba/your_reads.fastq.qz

Am rulat ex02_fastq_stats.py:
scriptul a fost rulat pentru accesion-ul descarcat SRR684066
am obtinut:
Reads: 5483467
Mean length: 46.00
N rate: 0.0002
Mean Phred: 31.37

pentru demo01_fastq_qc.py am folosit: CPCT12345678R_AHHKYHDSXX_S13_L001_R1_001.fastq.gz
pentru ex02_fastq_stats.py am folosit:  SRR684066

De ce este esențială verificarea calității datelor înainte de analiza variantelor?
1.Prevenim erorile false: daca fisierele FASTQ contin multe erori, acestea pot fi interpretate gresit ca variante genetice false, pot aparea false negatives din cauza calitatii slabe a bazelor.
2.Asigurarea acuratetei alinierii, datele gresite pot duce la aliniere incorecta a fragmentelor pe genomul de referinta.
3.Controlul calitatii permite detectarea dezechilibrarea intre probe, contaminarilor, distribuirii neuniforme a acoperirii genomului.
4.Analiza calitatitii furnizeaza informatii despre phred scores, lungimea fragmentelor, acoperirea, GC ce ghideaza alegerea parametrilor de filtrare in pasii urmatori
5.Asigura increderea in rezultate.
