import numpy as np
import pandas as pd
from pathlib import Path

handle = "indibotoc"
out_path = Path(f"data/work/{handle}/lab06/expression_matrix.csv")
out_path.parent.mkdir(parents=True, exist_ok=True)

N_GENES = 1000    
N_SAMPLES = 20    

# Gene IDs
genes = [f"Gene_{i+1}" for i in range(N_GENES)]
samples = [f"Sample_{j+1}" for j in range(N_SAMPLES)]

rng = np.random.default_rng(42)
expr = rng.lognormal(mean=4, sigma=1.0, size=(N_GENES, N_SAMPLES))
df = pd.DataFrame(expr, index=genes, columns=samples)

df.to_csv(out_path)
print(f"[OK] Am generat matricea de expresie: {out_path}")
print(df.head())