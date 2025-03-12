# Hadi Ahmad
# ahmad287
# 03/06/2025
# MyHPTest.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, ranksums
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

def my_hp_test(data, dataset_name):
    genes = list(data.keys())
    num_samples = len(next(iter(data.values()))) // 2
    group1 = {gene: values[:num_samples] for gene, values in data.items()}
    group2 = {gene: values[num_samples:] for gene, values in data.items()}
    
    # filter genes with 0 expression values across all the patients in the union of the two groups.
    filtered_genes = [gene for gene in genes if sum(data[gene]) > 0]
    
    p_values_ttest = []
    p_values_ranksum = []
    
    for gene in filtered_genes:
        _, p_t = ttest_ind(group1[gene], group2[gene], equal_var=False)
        _, p_r = ranksums(group1[gene], group2[gene])
        p_values_ttest.append((gene, p_t))
        p_values_ranksum.append((gene, p_r))
    
    # sort genes by p-value
    p_values_ttest.sort(key=lambda x: x[1])
    p_values_ranksum.sort(key=lambda x: x[1])
    
    # lists top 10 genes
    print(f"Top 10 genes by t-test for {dataset_name}:")
    for gene, p in p_values_ttest[:10]:
        print(f"{gene}: {p}")
    
    print(f"\nTop 10 genes by Wilcoxon rank-sum test for {dataset_name}:")
    for gene, p in p_values_ranksum[:10]:
        print(f"{gene}: {p}")
    
    # count significant genes
    significant_ttest = sum(1 for _, p in p_values_ttest if p < 0.05)
    significant_ranksum = sum(1 for _, p in p_values_ranksum if p < 0.05)
    print(f"\nNumber of significant genes in t-test: {significant_ttest}")
    print(f"Number of significant genes in Wilcoxon rank-sum test: {significant_ranksum}\n")
    
    # histogram plotted
    plt.figure(figsize=(12, 5))
    plt.hist([p for _, p in p_values_ttest], bins=50, alpha=0.5, label="t-test")
    plt.hist([p for _, p in p_values_ranksum], bins=50, alpha=0.5, label="Wilcoxon rank-sum")
    plt.xlabel("p-value")
    plt.ylabel("Frequency")
    plt.title(f"P-value Distribution for {dataset_name}")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python MyHPTest.py <input_file>")
        sys.exit(1)
    input_file = sys.argv[1]
    header, data = load_data(input_file)
    my_hp_test(data, input_file)
