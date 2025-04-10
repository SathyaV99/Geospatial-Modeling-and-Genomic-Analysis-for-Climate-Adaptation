import pandas as pd

def extract_genes_from_csv(csv_path, feature_types=["CDS", "gene", "mRNA"], save_as=None):
    """
    Extract gene-related entries from a .csv file converted from .gbff.

    Parameters:
    - csv_path: str — Path to the input CSV file.
    - feature_types: list — List of feature types to extract (default: CDS, gene, mRNA).
    - save_as: str — Optional output CSV filename to save results.

    Returns:
    - DataFrame with filtered gene-related entries.
    """
    df = pd.read_csv(csv_path)

    # Filter gene-related features
    gene_df = df[df['Feature_Type'].isin(feature_types)].copy()

    # Optional save
    if save_as:
        gene_df.to_csv(save_as, index=False)
        print(f"[✓] Extracted {len(gene_df)} {feature_types} features from {csv_path} to {save_as}")

    return gene_df

# === Paths to CSV files (Change these if your path differs) ===
base_path = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations - NONGITHUB\Gene_Feature_Extraction\1_genomic_feature_extraction"

# === Extracting CDS features for each species ===
wildyak_path = f"{base_path}\\wildyak_genomic_features.csv"
takin_path = f"{base_path}\\takin_genomic_features.csv"
waterbuffalo_path = f"{base_path}\\waterbuffalo_genomic_features.csv"

extract_genes_from_csv(wildyak_path, feature_types=["CDS"], save_as=f"{base_path}\\wildyak_CDS.csv")
extract_genes_from_csv(takin_path, feature_types=["CDS"], save_as=f"{base_path}\\takin_CDS.csv")
extract_genes_from_csv(waterbuffalo_path, feature_types=["CDS"], save_as=f"{base_path}\\waterbuffalo_CDS.csv")
