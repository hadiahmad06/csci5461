import sys
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.metrics import accuracy_score

def Q3_1_3(data_file):
    data = np.load(data_file)
    
    # extracts data
    X_train = data['training_data']
    y_train = data['training_label']
    X_test = data['testing_data']
    y_test = data['testing_label']
    
    # selects the top 1000 genes (first 1000 columns)
    X_train = X_train[:, :1000]
    X_test = X_test[:, :1000]
    
    # init and train linear SVM
    svm = LinearSVC(max_iter=10000)
    svm.fit(X_train, y_train)

    y_pred = svm.predict(X_test)
    
    # calculate accuract score
    accuracy = accuracy_score(y_test, y_pred)
    print(f"Classification accuracy: {accuracy:.4f}")
    
    return accuracy

if __name__ == "__main__":
    data_file = sys.argv[1]
    Q3_1_3(data_file)