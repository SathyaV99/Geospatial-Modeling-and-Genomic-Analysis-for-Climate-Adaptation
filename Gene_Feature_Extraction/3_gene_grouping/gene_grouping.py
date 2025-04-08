import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import csv

def detect_separator(file_path):
    """Auto-detect CSV delimiter (comma or tab)."""
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
    output_dir="gene_family_output",
    output_csv="gene_family_comparison.csv",
    output_png="gene_family_heatmap.png"
):
    """
    Analyze and compare gene family keyword counts across species.

    Args:
        csv_paths_dict (dict): {species_name: path_to_csv}
        keyword_list (list): List of keywords to search for in Product column
        output_dir (str): Directory to save output files
        output_csv (str): CSV filename for keyword comparison
        output_png (str): PNG filename for heatmap

    Returns:
        pd.DataFrame: Keyword count matrix (keywords x species)
    """
    if keyword_list is None:
        keywords = [
    # Original + GPCR family
    "kinase", "GPCR", "G protein-coupled", "7TM", "seven transmembrane", "zinc", "transport", "receptor","cytochrome", "ATPase", "ubiquitin", "oxidase", "synthetase",

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
    "growth factor", "morphogen", "developmental", "Wnt", "BMP", "Notch"]

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

            # Debug: Show available columns
            if species not in comparison_data:
                print(f"📋 Columns in {species}: {df.columns.tolist()}")

            if "feature_type" not in df.columns or "product" not in df.columns:
                print(f"⚠️ Missing required columns in {species} — skipping")
                continue

            cds_df = df[df["feature_type"] == "CDS"]

            counts = {}
            for keyword in keyword_list:
                count = cds_df["product"].str.contains(keyword, case=False, na=False).sum()
                counts[keyword] = count

            comparison_data[species] = counts

        except Exception as e:
            print(f"❌ Error processing {species}: {e}")

    result_df = pd.DataFrame(comparison_data).fillna(0).astype(int)

    if result_df.empty:
        print("❌ No data to visualize. Exiting.")
        return result_df

    # Save to CSV
    csv_path = os.path.join(output_dir, output_csv)
    result_df.to_csv(csv_path)
    print(f"✅ Saved keyword comparison to: {csv_path}")

    # Save very large heatmap
    plt.figure(figsize=(40, 20))  # Make the figure extremely large
    sns.heatmap(result_df, annot=True, cmap="YlGnBu", fmt="d", cbar_kws={'label': 'Count'})
    plt.title("Gene Family Keyword Counts Across Species", fontsize=24)
    plt.ylabel("Gene Family Keywords", fontsize=18)
    plt.xlabel("Species", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()

    png_path = os.path.join(output_dir, output_png)
    plt.savefig(png_path, dpi=300)  # High DPI for quality
    plt.close()
    print(f"✅ Saved large heatmap to: {png_path}")

    return result_df


######################################################################

csv_files = {
    "Takin": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\1_genomic_feature_extraction\takin_genomic_features.csv",
    "Wild Yak": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\1_genomic_feature_extraction\wildyak_genomic_features.csv",
    "Water Buffalo": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\1_genomic_feature_extraction\waterbuffalo_genomic_features.csv"
}


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


df_result = analyze_gene_families(csv_files, keyword_list=keywords)
print(df_result)
