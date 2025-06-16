import pandas as pd
from scipy.stats import kruskal
from glob import glob

# Cargar todos los CSV en una lista
csv_files = glob("errores*.csv")
dfs = [pd.read_csv(file) for file in csv_files]
df_all = pd.concat(dfs, ignore_index=True)

# Agrupar por parámetro
parametros = df_all["Parámetro"].unique()

print("RESULTADOS DEL TEST DE KRUSKAL-WALLIS POR PARÁMETRO:\n")

for param in parametros:
    subset = df_all[df_all["Parámetro"] == param]
    grupos = [grupo["Error"].values for _, grupo in subset.groupby("Algoritmo")]
    
    # Aplicar Kruskal-Wallis
    stat, p_value = kruskal(*grupos)
    print(f"Parámetro: {param}")
    print(f"  H = {stat:.4f}, p = {p_value:.4g}")
    if p_value < 0.05:
        print("Diferencias significativas entre algoritmos.\n")
    else:
        print("No hay diferencias significativas.\n")
