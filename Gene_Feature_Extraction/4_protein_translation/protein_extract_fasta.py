import pandas as pd

def extract_subset_fasta(csv_path, fasta_out, n=5000):
    df = pd.read_csv(csv_path, sep=None, engine='python')
    df.columns = [col.strip().lower() for col in df.columns]

    count = 0
    with open(fasta_out, "w") as f:
        for idx, row in df.iterrows():
            if pd.notna(row.get("protein_id")) and pd.notna(row.get("translation")):
                f.write(f">{row['protein_id']}\n{row['translation']}\n")
                count += 1
                if count >= n:
                    break

    print(f"✅ Saved {count} proteins to: {fasta_out}")


extract_subset_fasta(
    r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\1_genomic_feature_extraction\wildyak_genomic_features.csv",
    r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\4_protein_translation\protein_fastas\wild_yak_subset.fasta",
    n=3000
)
