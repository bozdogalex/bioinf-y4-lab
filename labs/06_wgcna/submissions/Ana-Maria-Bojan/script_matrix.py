import pandas as pd
import numpy as np
import os


NUM_GENES = 200       
NUM_SAMPLES = 6     
HANDLE = "Ana-Maria-Bojan"  
LAB_DIR = f"data/work/{HANDLE}/lab06/"
OUTPUT_FILENAME = "expression_matrix.csv"


real_genes = {
    # 1. GAPDH: Housekeeping (expresie mare si stabila)
    "GAPDH": np.random.uniform(low=1000.0, high=1200.0, size=NUM_SAMPLES),
    # 2. TP53: Tumor Suppressor (expresie medie)
    "TP53": np.random.uniform(low=50.0, high=80.0, size=NUM_SAMPLES),
    # 3. MYC: Proto-oncogene (expresie variabila, posibila supraexpresie)
    "MYC": np.concatenate([
        np.random.uniform(low=10.0, high=30.0, size=3), # 3 probe cu expresie mica
        np.random.uniform(low=150.0, high=250.0, size=3) # 3 probe cu supraexpresie simulata
    ])
}

# 1. Genereaza Genele Ramase ca Placeholder-uri
num_placeholder_genes = NUM_GENES - len(real_genes)
placeholder_gene_ids = [f'GENE_{i:04d}' for i in range(1, num_placeholder_genes + 1)]
all_gene_ids = list(real_genes.keys()) + placeholder_gene_ids

# 2. Genereaza Date de Expresie pentru Placeholder-uri
# Valori aleatorii simulate (folosind o distributie normala logaritmica)
placeholder_data = np.exp(np.random.normal(loc=4.0, scale=1.0, size=(num_placeholder_genes, NUM_SAMPLES)))

# 3. Combina Datele (Gene Reale + Placeholder-uri)
real_data_array = np.array(list(real_genes.values()))
combined_data = np.vstack([real_data_array, placeholder_data])
combined_data = np.round(combined_data, 2)

# 4. Creeaza DataFrame-ul
sample_ids = [f'GSM{i:03d}' for i in range(1, NUM_SAMPLES + 1)]
df = pd.DataFrame(combined_data, index=all_gene_ids, columns=sample_ids)

# 5. Salveaza in Format CSV (Indeplineste cerinta)
os.makedirs(LAB_DIR, exist_ok=True)
output_path = os.path.join(LAB_DIR, OUTPUT_FILENAME)
df.to_csv(output_path, index=True, header=True)

print("--- Matrice de Expresie de Test Generată ---")
print(f"Matricea a fost salvată la: {output_path}")
print(f"Gene incluse: {list(real_genes.keys())} + {num_placeholder_genes} gene placeholder")
print("\nPrimele 5 rânduri (inclusiv genele reale):")
print(df.head())