import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def compare_interpro_domains(tsv_files, output_dir="interpro_output", top_n=20):
    os.makedirs(output_dir, exist_ok=True)
    domain_counts = {}

    for species, path in tsv_files.items():
        print(f"📂 Processing: {species}")
        try:
            df = pd.read_csv(path, sep="\t", header=None)
            df.columns = [
                "Protein_ID", "MD5", "Length", "Analysis", "Signature_Acc",
                "Signature_Desc", "Start", "End", "Score", "Status", "Date",
                "InterPro_Acc", "InterPro_Desc", "GO", "Pathways"
            ][:len(df.columns)]

            top_domains = df["InterPro_Desc"].value_counts().head(top_n)
            domain_counts[species] = top_domains

        except Exception as e:
            print(f"❌ Error reading {species}: {e}")

    # Combine into one matrix
    if domain_counts:
        domain_df = pd.concat(domain_counts, axis=1).fillna(0).astype(int)
        domain_df.to_csv(os.path.join(output_dir, "interpro_domain_comparison.csv"))
        print(f"✅ Saved CSV: interpro_domain_comparison.csv")

        # Plot
        plt.figure(figsize=(15, 8))
        sns.heatmap(domain_df, annot=True, cmap="YlGnBu", fmt="d")
        plt.title("Top InterPro Domains Across Species")
        plt.ylabel("Domain Description")
        plt.xlabel("Species")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "interpro_domain_heatmap.png"))
        plt.close()
        print(f"✅ Saved heatmap: interpro_domain_heatmap.png")

# 🔧 Update with your real file paths
tsv_files = {
    "Water Buffalo": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\4_protein_translation\buffalo_interpro.tsv",
    "Wild Yak": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\4_protein_translation\wildyak_full_interpro.tsv",
    "Takin": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\4_protein_translation\takin_interpro.tsv"
}

compare_interpro_domains(tsv_files)
