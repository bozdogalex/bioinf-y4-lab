# Lab 07 – Network Viz Notes 

## Layout folosit

Pentru vizualizarea rețelei de co-expresie am folosit **layout-ul `spring` din NetworkX** (`nx.spring_layout`), cu un seed fix (SEED=42) pentru a obține o poziționare stabilă/reproducibilă a nodurilor.

Motivație:
- `spring_layout` este un algoritm de tip **force-directed**: nodurile care au multe muchii între ele tind să fie „apropiate” vizual, iar componentele/modulii se separă natural în spațiu.
- Fiind o rețea relativ mică (14 noduri după filtrări), layout-ul este rapid și scoate în evidență bine structura de module și poziția genelor hub.

## Reflecție: 

În Lab 6 am lucrat aproape exclusiv cu:
- **măsuri numerice** (corelații, praguri, număr de gene pe modul),
- și rezultate „textuale” (mapping gene → modul, număr de module detectate).

Vizualizarea rețelei în Lab 7 aduce câteva avantaje clare:

1. **Înțelegere structurală imediată**  
   Dintr-o singură imagine văd:
   - cum se grupează genele pe module (culori diferite),
   - cât de „dense” sunt modulele (câte muchii interne au),
   - dacă există gene care fac legătura între module (posibile „bridge genes”).

2. **Hub genes sunt mult mai evidente**  
   În Lab 6 vedeam doar un top numeric.  
   În rețea, hub-urile (ex. **DNASE1L3, FCN3, S100A9, C1QA/B/C**) apar:
   - cu **dimensiune mai mare a nodului** (proporțional cu degree),
   - poziționate central în modulul lor.  
   E mult mai ușor să apreciez vizual „rolul” acestor gene în topologia rețelei.

3. **Raport direct între module și conectivitate**  
   Numeric știu doar că există N module; vizual pot vedea:
   - care module sunt compacte și bine conectate,
   - care sunt mai dispersate,
   - dacă există module aproape complet izolate (potențial artefact sau grup funcțional foarte specific).

4. **Suport mai bun pentru interpretare biologică**  
   Când un modul mic, foarte dens, este dominat de gene cu rol imun/inflamator, imaginea ajută să argumentez că:
   - aceste gene formează un „subprogram” co-reglat,
   - hub-urile din modul ar putea fi ținte candidate pentru analize funcționale suplimentare.  

În concluzie, vizualizarea rețelei nu înlocuiește analiza numerică din Lab 6, dar:
- o **completează**,  
- permite verificarea intuitivă a structurii găsite (sanity check),  
- și face mult mai ușoară comunicarea rezultatelor (de exemplu într-un raport sau poster).
