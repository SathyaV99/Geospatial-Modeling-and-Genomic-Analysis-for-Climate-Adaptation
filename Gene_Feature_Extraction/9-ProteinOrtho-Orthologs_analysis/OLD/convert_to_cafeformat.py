import pandas as pd
import os

def convert_proteinortho_to_cafe(proteinortho_tsv, species_names, output_file):
    # Load proteinortho file
    df = pd.read_csv(proteinortho_tsv, sep='\t', header=None, comment='#')

    # Assign headers manually
    columns = ['Orthogroup', 'Genes', 'Species_with_gene'] + species_names
    df.columns = columns[:len(df.columns)]

    # Count number of genes per species (split by ",")
    gene_counts = df[species_names].fillna('*').applymap(lambda x: 0 if x == '*' else len(x.split(',')))

    # Add Orthogroup column
    gene_counts.insert(0, 'Orthogroup', df['Orthogroup'])

    # Save to file
    gene_counts.to_csv(output_file, sep='\t', index=False)
    print(f"[✓] CAFÉ input table saved to: {output_file}")

# === Run This ===
proteinortho_path = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\proteinortho\myproject.proteinortho.tsv"
output_path = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\cafe_input_gene_counts.tsv"

species = ['water_buffalo_proteins.fasta', 'wild_yak_proteins.fasta', 'takin_proteins.fasta']

convert_proteinortho_to_cafe(proteinortho_path, species, output_path)
