import pandas as pd
import os

def extract_core_orthologs_from_proteinortho(input_tsv, output_dir):
    """
    Extract core orthologs from a proteinortho .tsv file (no headers).
    Core = orthogroups present in all species (no missing entries).

    Parameters:
        input_tsv (str): Path to proteinortho .tsv file.
        output_dir (str): Output directory where CSV will be saved.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define headers manually (assuming 3 species)
    headers = ["Orthogroup_ID", "Connectivity", "Species_Count",
               "Water_Buffalo", "Wild_Yak", "Takin"]

    # Load the TSV with custom headers
    df = pd.read_csv(input_tsv, sep='\t', header=None, names=headers, comment="#")

    # Filter for rows where all 3 species have hits (i.e., core orthologs)
    core_df = df.dropna(subset=["Water_Buffalo", "Wild_Yak", "Takin"])

    # Save core orthologs
    output_path = os.path.join(output_dir, "core_orthologs_all_species.csv")
    core_df.to_csv(output_path, index=False)

    print(f"[✓] Extracted {len(core_df)} core orthologs.")
    print(f"[✓] Saved to: {output_path}")
    return core_df

# === Example usage ===
input_tsv_path = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\proteinortho\myproject.proteinortho.tsv"

output_folder = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\core_ortholog_output"

# Run it
extract_core_orthologs_from_proteinortho(input_tsv_path, output_folder)
