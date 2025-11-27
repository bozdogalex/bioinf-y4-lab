# Note Laborator 7 - Vizualizare Rețele și Diseasome
**Handle:** coroiusamy

## 1. Metoda de Layout Folosită
Am utilizat **Spring Layout** (`nx.spring_layout` din NetworkX).
- **Algoritm:** Fruchterman-Reingold (force-directed graph drawing).
- **Principiu:** Acest algoritm simulează un sistem fizic în care nodurile se comportă ca sarcini electrice (se resping reciproc), iar muchiile acționează ca arcuri (atrag nodurile conectate).
- **Rezultat:** Nodurile care sunt puternic interconectate (același modul) se grupează natural în clusteri vizuali, iar nodurile neconectate sunt împinse la distanță.

---

## 2. Reflecție: Vizualizare vs. Analiză Numerică (Lab 6)

În Laboratorul 6 am obținut fișiere CSV (liste de gene și module). Deși analiza numerică este precisă, vizualizarea din Laboratorul 7 aduce avantaje cruciale:

1.  **Percepția "Gestalt" (Imaginea de ansamblu):**
    - Analiza numerică ne spune că "Gena A și Gena B sunt în Modulul 1".
    - Vizualizarea ne arată **densitatea** conexiunilor. În imaginea generată, se vede clar cum nodurile verzi formează un cluster strâns (clica), sugerând o reglementare foarte strictă între ele (co-expresie puternică), în timp ce modulele sunt complet separate spațial.

2.  **Identificarea Topologiei Hub-urilor:**
    - Într-un tabel, un hub este doar un număr (Grad: 6).
    - Vizual, putem vedea **poziția strategică** a hub-ului. Vedem dacă un hub este "central" într-un cluster (cum este *APOA2* în centrul grupului verde) sau dacă este un "bridge" (punte) care leagă două module diferite (în cazul meu, modulele sunt disjuncte, ceea ce indică procese biologice independente).

3.  **Detectarea Erorilor sau Artefactelor:**
    - Vizualizarea a permis observarea imediată a faptului că rețeaua este formată din două componente conexe separate (una mare verde, una mică albastră/mov), lucru greu de dedus rapid doar citind o matrice de adiacență.

---

## 3. Bonus: Gephi

**Gephi:**
Am utilizat algoritmul **ForceAtlas2** pentru a spația nodurile. Diferența majoră a fost interactivitatea: am putut vedea în timp real cum cele două module se separă fizic din cauza lipsei de conexiuni între ele (insule distincte). Am ajustat parametrul "Gravity" pentru a aduce ambele componente pe același ecran, lucru mult mai intuitiv de făcut vizual decât prin cod.