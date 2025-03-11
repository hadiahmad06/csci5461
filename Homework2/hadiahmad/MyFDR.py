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

def my_fdr(filename, num_genes):
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
    
    # sort by p-values
    p_values_ttest.sort(key=lambda x: x[1])
    selected_genes = p_values_ttest[:num_genes]
    
    # extract p-values
    _, selected_p_values = zip(*selected_genes)
    
    # apply FDR correction
    _, fdr_corrected, _, _ = multipletests(selected_p_values, method='fdr_bh')
    
    # compute upper bound of FDR
    fdr_bound = max(fdr_corrected)
    
    print(f"Upper bound of FDR for selecting {num_genes} genes: {fdr_bound:.4f}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python MyFDR.py <input_file> <num_genes>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    num_genes = int(sys.argv[2])
    if num_genes not in [20, 50, 100, 200]:
        print("Error: Number of genes must be 20, 50, 100, or 200.")
        sys.exit(1)
    my_fdr(input_file, num_genes)