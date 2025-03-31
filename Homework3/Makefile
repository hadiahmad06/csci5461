PYTHON = python3

# data files
DATA1 = data/Data1.npz
DATA2 = data/Data2.npz

# scripts to run
Q1_SCRIPT = scripts/Q1.py
Q2_1_SCRIPT = scripts/Q2_1.py
Q3_1_2_SCRIPT = scripts/Q3_1_2.py
Q3_3_SCRIPT = scripts/Q3_3.py

all: 1 2 3

0:
	$(PYTHON) $(Q1_SCRIPT) $(DATA1)

1:
	@echo "Running K-means clustering with k=10, 20, 50..."
	$(PYTHON) $(Q2_1_SCRIPT) $(DATA1) 10
	$(PYTHON) $(Q2_1_SCRIPT) $(DATA1) 20
	$(PYTHON) $(Q2_1_SCRIPT) $(DATA1) 50

2:
	@echo "Running KNN classification with k=1, 3, 5 using all genes..."
	$(PYTHON) $(Q3_1_2_SCRIPT) $(DATA2) 1 1
	$(PYTHON) $(Q3_1_2_SCRIPT) $(DATA2) 3 1
	$(PYTHON) $(Q3_1_2_SCRIPT) $(DATA2) 5 1

	@echo "Running KNN classification with k=1, 3, 5 using top 1000 genes..."
	$(PYTHON) $(Q3_1_2_SCRIPT) $(DATA2) 1 0
	$(PYTHON) $(Q3_1_2_SCRIPT) $(DATA2) 3 0
	$(PYTHON) $(Q3_1_2_SCRIPT) $(DATA2) 5 0

3:
	@echo "Running SVM classification using top 1000 genes..."
	$(PYTHON) $(Q3_3_SCRIPT) $(DATA2)

clean:
	@echo "Cleaning up generated files..."
	rm -rf __pycache__ *.png *.log *.out

zip:
	@echo "Creating submission package..."
	zip -r hw3.zip report.pdf scripts/ data/

.PHONY: all test_q2 test_q3_knn test_q3_svm clean zip