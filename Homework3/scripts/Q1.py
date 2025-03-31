import sys
import numpy as np
import gseapy as gp

def main(data_file):
    # Load the .npz data file
    data = np.load(data_file)

    # Extract genes (update the key if needed)
    gene_key = 'Gene_Name'  # Change this if your file uses a different key
    if gene_key in data:
        genes = data[gene_key]
    else:
        print(f"Error: Key '{gene_key}' not found in the dataset.")
        return

    # Get the first 200 genes
    first_200_genes = genes[:200]

    # Find GO terms using gseapy
    print("\nRetrieving GO Terms...")
    enr = gp.enrichr(gene_list=list(first_200_genes), gene_sets='GO_Biological_Process_2021', organism='human')

    # Print the top 10 GO term codes
    if enr is not None and not enr.res2d.empty:
        print("\nTop 10 GO Term Codes:")
        for i, row in enr.res2d.head(10).iterrows():
            go_code = row['Term'].split('(')[-1].strip(')')  # Extract GO code from term
            print(go_code)
    else:
        print("\nNo GO terms found for the given genes.")

if __name__ == "__main__":
    
    main(sys.argv[1])