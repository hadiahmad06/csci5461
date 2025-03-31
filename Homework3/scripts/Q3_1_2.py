import sys
import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score

def Q3_1_2(data_file, k, flag):
    data = np.load(data_file)

    # extracts datasets
    gene_names = data['Gene_Name']
    X_train = data['training_data'].T
    y_train = data['training_label']
    X_test = data['testing_data'].T
    y_test = data['testing_label']

    # print(X_train.shape)
    # print(y_train.shape)

    # if flag is 0, use only the top 1000 genes
    if flag == 0:
        top_genes = gene_names[:1000]
        X_train = X_train[:1000, :]
        X_test = X_test[:1000, :]

    # transpose X_train and X_test to match (n_samples, n_features)
    X_train = X_train.T
    X_test = X_test.T

    # trains KNN model
    knn = KNeighborsClassifier(n_neighbors=k)
    knn.fit(X_train, y_train)

    # predict and compute accuracy
    y_pred = knn.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)

    print(f'Classification Accuracy (k={k}, flag={flag}): {accuracy:.4f}')
    return accuracy

if __name__ == "__main__":
    data_file = sys.argv[1]
    k = int(sys.argv[2])
    flag = int(sys.argv[3])
    Q3_1_2(data_file, k, flag)