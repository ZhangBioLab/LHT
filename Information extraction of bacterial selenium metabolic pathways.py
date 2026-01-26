import os
import pandas as pd


genomes_lht_dir = '/home/lihengtao/data/prokaryotes/original/genomes_new/genomes_lht'
all_bacteria_dir = '/home/lihengtao/data1/prokaryotes/all_bacteria_2025'


data = []


classifications = [d for d in os.listdir(all_bacteria_dir) if os.path.isdir(os.path.join(all_bacteria_dir, d))]


classification_stats = {}

for classification in classifications:
    classification_genomes_path = os.path.join(genomes_lht_dir, classification)
    classification_bacteria_path = os.path.join(all_bacteria_dir, classification)
    
    
    if not os.path.isdir(classification_genomes_path):
        continue

    
    accessions = [d for d in os.listdir(classification_genomes_path) if os.path.isdir(os.path.join(classification_genomes_path, d))]
    total_accessions = len(accessions) 

    
    final_dir = os.path.join(classification_bacteria_path, 'final')
    if not os.path.exists(final_dir):
        continue
    proteins = [d for d in os.listdir(final_dir) if os.path.isdir(os.path.join(final_dir, d))]
    
    proteins = sorted(proteins)

    
    protein_counts = {protein: 0 for protein in proteins}

    for accession in accessions:
        accession_data = {'classifications': classification, 'accessions': accession}
        for protein in proteins:
            fasta_file = os.path.join(final_dir, protein, accession + '.fasta')
            if os.path.exists(fasta_file):
                
                count = 0
                with open(fasta_file, 'r') as f:
                    for line in f:
                        if line.startswith('>'):
                            count += 1
                if count >= 2:
                    accession_data[protein] = f'+({count})'
                else:
                    accession_data[protein] = '+'
                
                protein_counts[protein] += 1
            else:
                accession_data[protein] = '-'
        data.append(accession_data)
    
    
    classification_name_with_count = f"{classification}({total_accessions})"
    classification_stats[classification_name_with_count] = protein_counts


df = pd.DataFrame(data)


columns = ['classifications', 'accessions'] + proteins
df = df[columns]


stats_data = []
for classification_name, protein_counts in classification_stats.items():
    row = {'classifications(count)': classification_name}
    row.update(protein_counts)
    stats_data.append(row)

stats_df = pd.DataFrame(stats_data)


stats_columns = ['classifications(count)'] + proteins
stats_df = stats_df[stats_columns]


with pd.ExcelWriter('./statistical_results_2025.xlsx') as writer:
    df.to_excel(writer, sheet_name='all_accessions', index=False)
    stats_df.to_excel(writer, sheet_name='classification_statistics', index=False)

