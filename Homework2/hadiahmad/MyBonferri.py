# Hadi Ahmad
# ahmad287
# 03/06/2025
# MyBonferri.py

import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import sys

def load_data(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    header = lines[0].strip().split('\t')[1:]
    data = {}

    for line in lines[1:]:
        parts = line.strip().split('\t')
        gene = parts[0]
        values = list(map(float, parts[1:]))
        data[gene] = values
    
    return header, data

def my_bonferroni(filename):
    _, data = load_data(filename)
    genes = list(data.keys())
    num_samples = len(next(iter(data.values()))) // 2

    if num_samples < 2:
        print("Error: Not enough samples to perform t-test.", file=sys.stderr)
        return

    # remove genes with all zero expression
    filtered_genes = [gene for gene in genes if any(data[gene])]

    gene_p_values = []

    for gene in filtered_genes:
        group1 = data[gene][:num_samples]
        group2 = data[gene][num_samples:]

        # check for zero variance in either group
        if np.var(group1) == 0 or np.var(group2) == 0:
            continue  

        # check for nan value
        _, p_value = ttest_ind(group1, group2, equal_var=False)
        if np.isnan(p_value):  
            continue
        gene_p_values.append((gene, p_value))

    if not gene_p_values:
        print("No valid genes available for statistical testing.")
        return

    # sort before unpacking
    gene_p_values.sort(key=lambda x: x[1])
    genes_sorted, p_values = zip(*gene_p_values)

    # multiple testing correction (Bonferroni)
    _, bonferroni_corrected, _, _ = multipletests(p_values, method='bonferroni')

    # select significant genes
    selected_genes = [gene for gene, p in zip(genes_sorted, bonferroni_corrected) if p < 0.05]

    print(len(selected_genes))
    if selected_genes:
        print("\n".join(selected_genes))
    else:
        print("No significant genes found after Bonferroni correction.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python MyBonferroni.py <input_file>")
        sys.exit(1)

    my_bonferroni(sys.argv[1])