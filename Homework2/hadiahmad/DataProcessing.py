# # Hadi Ahmad
# # ahmad287
# # 03/06/2025
# # DataProcessing.py

import csv
import os

print(os.getcwd())
# Open the files
with open('./data/HiSeqV2.tsv', 'r') as rna_file, open('./data/HT_HG-U133A.tsv', 'r') as array_file, open('./data/ov_tcga_clinical_data.tsv', 'r') as clinical_file:
    
    # from RNAseq data
    rna_reader = csv.reader(rna_file, delimiter='\t')
    rna_header = next(rna_reader)
    rna_data = {header: [] for header in rna_header[1:]}
    # transposes data
    for row in rna_reader:
        for i, value in enumerate(row[1:], start=1):
            rna_data[rna_header[i]].append(value)

    # from microarray data
    array_reader = csv.reader(array_file, delimiter='\t')
    array_header = next(array_reader)
    # transposes data
    array_data = {header: [] for header in array_header[1:]}
    for row in array_reader:
        for i, value in enumerate(row[1:], start=1):
            array_data[array_header[i]].append(value)

    # from clinical data
    clinical_reader = csv.reader(clinical_file, delimiter='\t')
    clinical_header = next(clinical_reader)

    # gets indices of columns that will be checked
    patient_id_index = clinical_header.index('Sample ID')
    survival_status_index = clinical_header.index('Overall Survival Status')
    overall_survival_index = clinical_header.index('Overall Survival (Months)')

    rna_seq_output = []
    array_data_output = []

    for row in clinical_reader:
        sample_id = row[patient_id_index]
        survival_status = row[survival_status_index]
        overall_survival_str = row[overall_survival_index]

        # skip if non-float value
        if overall_survival_str == "NA":
            continue
        try:
            overall_survival = float(overall_survival_str)
        except ValueError:
            continue

        if survival_status == '1:DECEASED' and overall_survival < 36:
            group_id = '1'  # Group 1: DECEASED and survival < 36 months
        elif overall_survival > 36:
            group_id = '2'  # Group 2: survival > 36 months
        else:
            continue

        if sample_id in rna_header and sample_id in array_header:
            rna_row = [sample_id, group_id] + rna_data[sample_id]
            array_row = [sample_id, group_id] + array_data[sample_id]

            rna_seq_output.append(rna_row)
            array_data_output.append(array_row)

    with open('./output/SeqData.txt', 'w', newline='') as seq_file:
        seq_writer = csv.writer(seq_file, delimiter='\t')
        seq_writer.writerow(['sample_id', 'group_id'] + rna_header[1:])
        seq_writer.writerows(rna_seq_output)

    with open('./output/ArrayData.txt', 'w', newline='') as array_file:
        array_writer = csv.writer(array_file, delimiter='\t')
        array_writer.writerow(['sample_id', 'group_id'] + array_header[1:])
        array_writer.writerows(array_data_output)

    print("Preprocessing completed and data has been saved to 'SeqData.txt' and 'ArrayData.txt'")