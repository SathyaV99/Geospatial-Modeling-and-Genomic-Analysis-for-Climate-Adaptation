import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pprint
import os

def analyze_gbff_csv(csv_path, species_name):
    df = pd.read_csv(csv_path)
    summary = {}

    print(f"\n📊 Analysis for {species_name}")
    print("=" * 50)

    # Feature Type Counts
    feature_counts = df['Feature_Type'].value_counts()
    print("\n🔍 Feature Type Counts:")
    print(feature_counts)
    summary["feature_type_counts"] = feature_counts.to_dict()
    summary["unique_feature_types"] = list(feature_counts.index)

    # Gene & Protein Stats
    summary["total_rows"] = len(df)
    summary["columns"] = df.columns.tolist()
    summary["unique_genes"] = df['Gene'].nunique()
    summary["unique_proteins"] = df['Protein_ID'].nunique()
    print("\n🧬 Unique Genes:", summary["unique_genes"])
    print("🧬 Unique Proteins:", summary["unique_proteins"])

    # Locus Tags
    summary["unique_locus_tags"] = df['Locus_Tag'].nunique()
    print("\n🏷️ Unique Locus Tags:", summary["unique_locus_tags"])

    # Contig Stats
    summary["unique_contigs"] = df['Contig'].nunique()
    summary["contig_names"] = df['Contig'].unique().tolist()[:10]
    print("\n📦 Unique Contigs:", summary["unique_contigs"])

    # Strand Distribution
    summary["strand_distribution"] = df['Strand'].value_counts().to_dict()

    # CDS Length Stats
    cds = df[df['Feature_Type'] == 'CDS'].copy()
    cds['Length'] = cds['End'] - cds['Start']
    cds_lengths = cds['Length'].describe()
    summary["avg_gene_length"] = cds_lengths['mean']
    summary["longest_gene"] = cds_lengths['max']
    summary["shortest_gene"] = cds_lengths['min']
    print("\n📏 CDS Length Summary:")
    print(cds_lengths)

    # Plot CDS Length Distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(cds['Length'], bins=100, kde=True)
    plt.title(f'{species_name} - CDS Length Distribution')
    plt.xlabel('Length (bp)')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plot_path = os.path.splitext(csv_path)[0] + "_cds_plot.png"
    plt.savefig(plot_path)
    plt.close()
    print(f"📈 Saved CDS plot to: {plot_path}")

    # Top Products
    top_products = df['Product'].value_counts().head(10)
    summary["gene_product_counts"] = top_products.to_dict()
    print("\n🔥 Top Products by Count:")
    print(top_products)

    # Top Genes
    top_genes = df['Gene'].value_counts().head(10)
    summary["top_genes"] = top_genes.to_dict()
    print("\n🧪 Most Frequent Genes:")
    print(top_genes)

    # Missing product field count
    summary["missing_product_count"] = df['Product'].isnull().sum()

    # Translation lengths
    summary["top_translation_lengths"] = df['Translation'].dropna().apply(len).sort_values(ascending=False).head(5).tolist()

    # Notes Summary
    summary["notes_summary"] = df['Note'].dropna().value_counts().head(5).to_dict()

    # Non-coding features: tRNA, rRNA, misc_RNA
    for elem in ['tRNA', 'rRNA', 'misc_RNA']:
        subset = df[df['Feature_Type'] == elem]
        summary[f"{elem}_count"] = len(subset)
        if not subset.empty:
            print(f"\n🔬 {elem} count: {len(subset)}")
            print(subset[['Contig', 'Start', 'End', 'Strand', 'Product']].head())

    # Save summary to TXT
    summary_path = os.path.splitext(csv_path)[0] + "_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"📊 Summary Report for {species_name}\n")
        f.write("=" * 60 + "\n")
        f.write(pprint.pformat(summary, sort_dicts=False, width=120))
    print(f"\n📝 Saved summary report to: {summary_path}")
    print("=" * 50)
    print("✅ Analysis Completed")


# ▶️ Run for all
analyze_gbff_csv(
    r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\takin_genomic_features.csv",
    "Takin"
)

analyze_gbff_csv(
    r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\waterbuffalo_genomic_features.csv",
    "WaterBuffalo"
)

analyze_gbff_csv(
    r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations\Gene_Feature_Extraction\wildyak_genomic_features.csv",
    "WildYak"
)

