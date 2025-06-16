import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
from itertools import combinations

# Cargar todos los CSV 
from glob import glob
csv_files = glob("errores*.csv")
dfs = [pd.read_csv(file) for file in csv_files]
df_all = pd.concat(dfs, ignore_index=True)

# Obtener parámetros únicos
parametros = df_all["Parámetro"].unique()

print("RESULTADOS DEL TEST DE WILCOXON:\n")

for param in parametros:
    subset = df_all[df_all["Parámetro"] == param]
    pivot = subset.pivot(index="Ejecución", columns="Algoritmo", values="Error")
    
    algoritmos = pivot.columns.tolist()
    pares = list(combinations(algoritmos, 2))
    
    p_vals = []
    comparaciones = []

    for a1, a2 in pares:
        # Wilcoxon para muestras pareadas
        try:
            stat, p = wilcoxon(pivot[a1], pivot[a2])
            p_vals.append(p)
            comparaciones.append((a1, a2))
        except ValueError:
            # Si hay errores por muestras vacías o iguales, ignorar
            continue

    # Corrección Holm-Bonferroni
    rechazos, pvals_corr, _, _ = multipletests(p_vals, alpha=0.05, method='holm')

    print(f"\nParámetro: {param}")
    for i, ((a1, a2), p, p_corr, r) in enumerate(zip(comparaciones, p_vals, pvals_corr, rechazos)):
        signif = "si" if r else "no"
        print(f"  {a1} vs {a2}: p = {p:.4g}, p_corr = {p_corr:.4g} es {'diferencia significativa' if r else 'no significativa'} {signif}")