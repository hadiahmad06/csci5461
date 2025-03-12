import matplotlib.pyplot as plt

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

# function to get top selected genes based on average value
def get_top_genes(data, num_genes):
    # sort genes by the average value
    sorted_genes = sorted(data.items(), key=lambda x: sum(x[1])/len(x[1]), reverse=True)
    top_genes = {gene: values for gene, values in sorted_genes[:num_genes]}
    return top_genes

# function to get common genes between RNAseq and microarray data
def get_common_genes(rna_data, array_data):
    common_genes = set(rna_data.keys()).intersection(set(array_data.keys()))
    return len(common_genes)

if __name__ == "__main__":
    rna_filename = './data/HiSeqV2.tsv'
    array_filename = './data/HT_HG-U133A.tsv'
    
    header_rna, rna_data = load_data(rna_filename)
    header_array, array_data = load_data(array_filename)
    
    x_vals = []  # number of selected genes
    y_vals = []  # number of common genes

    for num_genes in range(10, 1001, 10):
        rna_top_genes = get_top_genes(rna_data, num_genes)
        array_top_genes = get_top_genes(array_data, num_genes)
        
        common_genes_count = get_common_genes(rna_top_genes, array_top_genes)
        
        x_vals.append(num_genes)
        y_vals.append(common_genes_count)
    
    # plot the data
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, y_vals, marker='.', linestyle='-', color='black')
    plt.title('Overlap of Gene Expression in RNAseq and Microarray Data')
    plt.xlabel('Number of Selected Genes')
    plt.ylabel('Number of Common Genes')
    plt.grid(True)
    plt.xticks(range(0, 1001, 100))
    plt.show()