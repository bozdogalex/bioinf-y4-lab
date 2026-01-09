- ce metodă de layout ați folosit (ex: spring, kamada-kawai),
    - Metoda de layout folosită a fost **Spring**.
- o scurtă reflecție: **Ce avantaje aduce vizualizarea față de analiza numerică din Lab 6?** 
    - În timp ce analiza numerică oferă rigoare statistică (p-values, coeficienți de corelație​), vizualizarea aduce o perspectivă vitală asupra arhitecturii sistemului:

        - **Percepția topologiei globale**: Tabelele cu mii de gene sunt abstracte. Vizualizarea grafică (ex. network plots, heatmaps) permite identificarea instantanee a structurii modulare („clusterelor”) și a gradului de separare dintre procesele biologice, lucruri greu de dedus doar din matricea de adiacență.

        - **Contextualizarea rolului genelor (Hub vs. Connector)**: Numeric, o genă poate avea un grad mare de conectivitate, dar vizualizarea ne arată unde este plasată: hub central sau punte care conectează două module diferite. Această poziționare spațială sugerează funcția biologică de integrare a semnalelor.

        - **Detectarea rapidă a erorilor**: Vizualizările sunt instrumente indispensabile tabelelor pentru Quality Control (QC). Ele permit observarea imediată a efectelor de batch, a outlier-ilor sau a artefactelor de clustering care ar putea distorsiona rezultatele numerice.