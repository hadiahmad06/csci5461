# Hadi Ahmad
# ahmad287
# 03/06/2025
# DataProcessing.py

import csv
import os

# print(os.getcwd())

with open('./data/HiSeqV2.tsv', 'r') as rna_file, open('./data/HT_HG-U133A.tsv', 'r') as array_file, open('./data/ov_tcga_clinical_data.tsv', 'r') as clinical_file:
    
    # from rnaseq data
    rna_reader = csv.reader(rna_file, delimiter='\t')
    rna_header = next(rna_reader)
    samples = rna_header[1:]  # Store patient/sample IDs
    genes = []  # Store gene names
    rna_data = {}

    for row in rna_reader:
        gene = row[0]
        genes.append(gene)
        rna_data[gene] = row[1:]

    # from microarray data
    array_reader = csv.reader(array_file, delimiter='\t')
    array_header = next(array_reader)
    array_data = {}

    for row in array_reader:
        gene = row[0]
        array_data[gene] = row[1:]

    clinical_reader = csv.reader(clinical_file, delimiter='\t')
    clinical_header = next(clinical_reader)

    # Get indices of relevant columns
    patient_id_index = clinical_header.index('Sample ID')
    survival_status_index = clinical_header.index('Overall Survival Status')
    overall_survival_index = clinical_header.index('Overall Survival (Months)')

    selected_samples = {}

    for row in clinical_reader:
        sample_id = row[patient_id_index]
        survival_status = row[survival_status_index]
        overall_survival_str = row[overall_survival_index]

        if overall_survival_str == "NA":
            continue
        try:
            overall_survival = float(overall_survival_str)
        except ValueError:
            continue

        if survival_status == '1:DECEASED' and overall_survival < 36:
            selected_samples[sample_id] = '1'  # Group 1: DECEASED and survival < 36 months
        elif overall_survival > 36:
            selected_samples[sample_id] = '2'  # Group 2: survival > 36 months

    rna_seq_output = [['sample_id', 'group_id'] + genes]
    for sample_id, group_id in selected_samples.items():
        if sample_id in samples:
            rna_seq_output.append([sample_id, group_id] + [rna_data[gene][samples.index(sample_id)] for gene in genes])

    array_output = [['sample_id', 'group_id'] + list(array_data.keys())]
    for sample_id, group_id in selected_samples.items():
        if sample_id in samples:
            array_output.append([sample_id, group_id] + [array_data[gene][samples.index(sample_id)] for gene in array_data])

    with open('./output/SeqData.txt', 'w', newline='') as seq_file:
        seq_writer = csv.writer(seq_file, delimiter='\t')
        seq_writer.writerows(rna_seq_output)

    with open('./output/ArrayData.txt', 'w', newline='') as array_file:
        array_writer = csv.writer(array_file, delimiter='\t')
        array_writer.writerows(array_output)

    print("Preprocessing completed and data has been saved to 'SeqData.txt' and 'ArrayData.txt'")