import pandas as pd
import os

def extract_protein_fastas(csv_paths_dict, output_dir="protein_fastas"):
    """
    Extracts protein FASTA files from .csv annotation files for multiple species.

    Args:
        csv_paths_dict (dict): Dictionary of {species_name: csv_path}
        output_dir (str): Output directory to save FASTA files
    """
    os.makedirs(output_dir, exist_ok=True)

    for species, csv_path in csv_paths_dict.items():
        try:
            df = pd.read_csv(csv_path, sep=None, engine='python')  # Auto-detect comma/tab
            df.columns = [col.strip().lower() for col in df.columns]

            if "protein_id" not in df.columns or "translation" not in df.columns:
                print(f"⚠️ Missing required columns in {species}. Skipping.")
                continue

            fasta_file = os.path.join(output_dir, f"{species.lower().replace(' ', '_')}_proteins.fasta")

            with open(fasta_file, "w") as f:
                count = 0
                for idx, row in df.iterrows():
                    protein_id = row["protein_id"]
                    sequence = row["translation"]
                    if pd.notna(protein_id) and pd.notna(sequence):
                        f.write(f">{protein_id}\n{sequence}\n")
                        count += 1

            print(f"✅ Saved {count} protein sequences for {species} → {fasta_file}")

        except Exception as e:
            print(f"❌ Error processing {species}: {e}")

csv_files = {
    "Takin": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\1_genomic_feature_extraction\takin_genomic_features.csv",
    "Wild Yak": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\1_genomic_feature_extraction\wildyak_genomic_features.csv",
    "Water Buffalo": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\1_genomic_feature_extraction\waterbuffalo_genomic_features.csv"
}

extract_protein_fastas(csv_files)

