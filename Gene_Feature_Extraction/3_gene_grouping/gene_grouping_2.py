import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import csv

def detect_separator(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        sample = f.read(1024)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=[',', '\t'])
            return dialect.delimiter
        except csv.Error:
            print(f"⚠️ Could not detect delimiter for {file_path}, defaulting to comma.")
            return ','

def analyze_gene_families(
    csv_paths_dict,
    keyword_list=None,
    keyword_groups=None,
    output_dir="gene_family_output",
    output_csv="gene_family_comparison.csv",
    grouped_csv="gene_function_summary.csv",
    output_png="gene_family_heatmap.png",
    grouped_png="gene_function_grouped_heatmap.png"
):
    if keyword_list is None:
        print("❌ No keywords provided.")
        return

    os.makedirs(output_dir, exist_ok=True)
    comparison_data = {}

    for species, path in csv_paths_dict.items():
        if not os.path.exists(path):
            print(f"⚠️ File not found: {path} — skipping {species}")
            continue

        try:
            sep = detect_separator(path)
            df = pd.read_csv(path, sep=sep, encoding="utf-8")
            df.columns = [col.strip().lower() for col in df.columns]

            if "feature_type" not in df.columns or "product" not in df.columns:
                print(f"⚠️ Missing required columns in {species} — skipping")
                continue

            cds_df = df[df["feature_type"] == "CDS"]
            counts = {kw: cds_df["product"].str.contains(kw, case=False, na=False).sum() for kw in keyword_list}
            comparison_data[species] = counts

        except Exception as e:
            print(f"X Error processing {species}: {e}")

    result_df = pd.DataFrame(comparison_data).fillna(0).astype(int)

    # Save keyword-level CSV and heatmap
    keyword_csv = os.path.join(output_dir, output_csv)
    result_df.to_csv(keyword_csv)
    print(f":D Saved keyword comparison to: {keyword_csv}")

    plt.figure(figsize=(40, 20))
    sns.heatmap(result_df, annot=True, cmap="YlGnBu", fmt="d", cbar_kws={'label': 'Count'})
    plt.title("Gene Family Keyword Counts Across Species", fontsize=24)
    plt.ylabel("Gene Family Keywords", fontsize=18)
    plt.xlabel("Species", fontsize=18)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_png), dpi=300)
    plt.close()
    print(f":D Saved keyword-level heatmap to: {output_png}")

    # Grouped analysis by shared function
    if keyword_groups:
        grouped_data = {}
        for group_name, keywords in keyword_groups.items():
            relevant = result_df.loc[result_df.index.intersection(keywords)]
            grouped_data[group_name] = relevant.sum()

        grouped_df = pd.DataFrame(grouped_data).T.astype(int)

        # Save grouped CSV and heatmap
        grouped_csv_path = os.path.join(output_dir, grouped_csv)
        grouped_df.to_csv(grouped_csv_path)
        print(f":D Saved grouped function summary to: {grouped_csv_path}")

        plt.figure(figsize=(12, 8))
        sns.heatmap(grouped_df, annot=True, cmap="YlGnBu", fmt="d", cbar_kws={'label': 'Count'})
        plt.title("Gene Function Category Counts Across Species", fontsize=20)
        plt.ylabel("Shared Function", fontsize=16)
        plt.xlabel("Species", fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, grouped_png), dpi=300)
        plt.close()
        print(f":D Saved grouped function heatmap to: {grouped_png}")

    return result_df

# Keyword list
keywords = [
    # Original + GPCR family
    "kinase", "GPCR", "G protein-coupled", "7TM", "seven transmembrane", "zinc", "transport", "receptor",
    "cytochrome", "ATPase", "ubiquitin", "oxidase", "synthetase",
    # Transcription / Regulation
    "transcription", "homeobox", "helix-loop-helix", "forkhead", "nuclear receptor", "TF", "coactivator", "corepressor",
    # Signal Transduction
    "phosphatase", "second messenger", "MAPK", "cAMP", "PKA", "calmodulin",
    # Stress / Detox
    "reductase", "peroxidase", "heat shock", "HSP", "proteasome", "stress",
    # Energy Metabolism
    "dehydrogenase", "glycolysis", "mitochondrial",
    # Hypoxia / Altitude
    "HIF", "hypoxia", "VEGF", "carbonic anhydrase", "oxygen", "myoglobin", "hemoglobin",
    # Immune Genes
    "interleukin", "toll-like", "TLR", "major histocompatibility", "MHC", "defensin", "immunoglobulin", "cytokine",
    # Development
    "growth factor", "morphogen", "developmental", "Wnt", "BMP", "Notch"
]

# Group keywords into biological functions
keyword_groups = {
    "Energy metabolism": ["ATPase", "cytochrome", "dehydrogenase", "mitochondrial"],
    "Protein quality control": ["ubiquitin", "proteasome", "heat shock", "oxidase"],
    "Signal transduction": ["kinase", "receptor", "calmodulin"],
    "Transcriptional regulation": ["zinc", "TF", "homeobox", "coactivator", "corepressor"],
    "Immune response": ["cytokine", "interleukin", "TLR", "MHC"],
    "Growth & development": ["Notch", "Wnt", "BMP", "growth factor"]
}

# CSV input
csv_files = {
    "Takin": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\1_genomic_feature_extraction\takin_genomic_features.csv",
    "Wild Yak": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\1_genomic_feature_extraction\wildyak_genomic_features.csv",
    "Water Buffalo": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\1_genomic_feature_extraction\waterbuffalo_genomic_features.csv"
}

# Run the analysis
analyze_gene_families(
    csv_paths_dict=csv_files,
    keyword_list=keywords,
    keyword_groups=keyword_groups
)
