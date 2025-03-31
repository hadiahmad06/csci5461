import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

def load_data(filename):
    # extracts the first 1000 gene expressions
    data = np.load(filename)
    gene_names = data["Gene_Name"][:1000]
    seq_data = data["SeqData"][:1000, :]

    return gene_names, seq_data

def kmeans_clustering(seq_data, k):
    # does k-means clustering and returns cluster labels
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(seq_data)
    return labels

def plot_histogram(labels, k):
    # plot histogram of clusters
    plt.hist(labels, bins=np.arange(-0.5, k + 0.5, 1), edgecolor="black")
    plt.xlabel("Cluster ID")
    plt.ylabel("Number of Genes")
    plt.title(f"K-means Clustering (k={k})")
    plt.xticks(range(k))
    plt.show()

def extract_cluster_genes(gene_names, labels, cluster_id):
    # extract genes belonging to a specific cluster
    selected_genes = gene_names[labels == cluster_id]
    print(f"Genes in Cluster {cluster_id}:")
    for gene in selected_genes:
        print(gene)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python Q2_1.py <data_file> k")
        sys.exit(1)

    data_file = sys.argv[1]
    k = int(sys.argv[2])

    gene_names, seq_data = load_data(data_file)
    labels = kmeans_clustering(seq_data, k)
    plot_histogram(labels, k)

    if k == 20:
        cluster_id = 0
        extract_cluster_genes(gene_names, labels, cluster_id)