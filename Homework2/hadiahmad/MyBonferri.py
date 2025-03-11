import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import sys

# load data
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
    group1 = {gene: values[:num_samples] for gene, values in data.items()}
    group2 = {gene: values[num_samples:] for gene, values in data.items()}
    
    # remove genes with zero expression across all patients
    filtered_genes = [gene for gene in genes if sum(data[gene]) > 0]
    
    p_values_ttest = []
    
    for gene in filtered_genes:
        _, p_t = ttest_ind(group1[gene], group2[gene], equal_var=False)
        p_values_ttest.append((gene, p_t))
    
    # multiple testing correction (Bonferroni)
    genes_sorted, p_values = zip(*p_values_ttest)
    _, bonferroni_corrected, _, _ = multipletests(p_values, method='bonferroni')
    
    # select significant genes
    selected_genes = [gene for gene, p in zip(genes_sorted, bonferroni_corrected) if p < 0.05]

    print(len(selected_genes))
    for gene in selected_genes:
        print(gene)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python MyBonferri.py <input_file> <num_genes>")
        sys.exit(1)
    my_bonferroni(sys.argv[1])
