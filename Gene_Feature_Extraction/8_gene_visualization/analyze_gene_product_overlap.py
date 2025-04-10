import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import seaborn as sns
import os

def analyze_gene_product_overlap(wildyak_csv, takin_csv, buffalo_csv, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Load CSVs
    wildyak_df = pd.read_csv(wildyak_csv)
    takin_df = pd.read_csv(takin_csv)
    buffalo_df = pd.read_csv(buffalo_csv)

    # Drop NAs and convert to sets
    wildyak_set = set(wildyak_df['Product'].dropna())
    takin_set = set(takin_df['Product'].dropna())
    buffalo_set = set(buffalo_df['Product'].dropna())

    # Shared and unique sets
    shared_all = wildyak_set & takin_set & buffalo_set
    unique_yak = wildyak_set - (takin_set | buffalo_set)
    unique_takin = takin_set - (wildyak_set | buffalo_set)
    unique_buffalo = buffalo_set - (wildyak_set | takin_set)
    only_yak = wildyak_set - buffalo_set
    only_buffalo = buffalo_set - wildyak_set
    only_takin = takin_set - buffalo_set

    # Save each gene product set to CSV
    pd.Series(list(shared_all)).to_csv(os.path.join(output_dir, "shared_all_species.csv"), index=False, header=["Gene_Product"])
    pd.Series(list(unique_yak)).to_csv(os.path.join(output_dir, "unique_to_wildyak.csv"), index=False, header=["Gene_Product"])
    pd.Series(list(unique_takin)).to_csv(os.path.join(output_dir, "unique_to_takin.csv"), index=False, header=["Gene_Product"])
    pd.Series(list(unique_buffalo)).to_csv(os.path.join(output_dir, "unique_to_waterbuffalo.csv"), index=False, header=["Gene_Product"])
    pd.Series(list(only_yak)).to_csv(os.path.join(output_dir, "only_yak_not_buffalo.csv"), index=False, header=["Gene_Product"])
    pd.Series(list(only_buffalo)).to_csv(os.path.join(output_dir, "only_buffalo_not_yak.csv"), index=False, header=["Gene_Product"])
    pd.Series(list(only_takin)).to_csv(os.path.join(output_dir, "only_takin_not_buffalo.csv"), index=False, header=["Gene_Product"])

    # Summary table
    summary = {
        "Shared (All 3)": len(shared_all),
        "Unique to Wild Yak": len(unique_yak),
        "Unique to Takin": len(unique_takin),
        "Unique to Water Buffalo": len(unique_buffalo),
        "Only in Wild Yak (not Buffalo)": len(only_yak),
        "Only in Water Buffalo (not Yak)": len(only_buffalo),
        "Only in Takin (not Buffalo)": len(only_takin),
    }

    # --- Save combined gene product presence matrix ---
    all_products = list(wildyak_set | takin_set | buffalo_set)
    matrix = pd.DataFrame({
        'Gene_Product': all_products,
        'In_WildYak': [p in wildyak_set for p in all_products],
        'In_Takin': [p in takin_set for p in all_products],
        'In_WaterBuffalo': [p in buffalo_set for p in all_products],
    })
    matrix.to_csv(os.path.join(output_dir, "gene_product_species_matrix.csv"), index=False)

    # --- Save shared gene products ---
    shared_genes_df = matrix[
        (matrix["In_WildYak"]) & (matrix["In_Takin"]) & (matrix["In_WaterBuffalo"])
    ]
    shared_genes_df.to_csv(os.path.join(output_dir, "gene_product_shared_all_species.csv"), index=False)

    # --- Venn Diagram ---
    plt.figure(figsize=(8, 6))
    venn3([wildyak_set, takin_set, buffalo_set], set_labels=('Wild Yak', 'Takin', 'Water Buffalo'))
    plt.title("Gene Product Overlap Across Species")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gene_product_venn.png"))
    plt.close()

    # --- Bar Chart ---
    summary_df = pd.DataFrame(list(summary.items()), columns=["Category", "Count"])
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Category", y="Count", data=summary_df, palette="viridis")
    plt.xticks(rotation=45, ha="right")
    plt.title("Gene Product Distribution Across Species")
    plt.ylabel("Number of Gene Products")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gene_product_distribution_barplot.png"))
    plt.close()

    # --- Heatmap ---
    heat_data = matrix[['In_WildYak', 'In_Takin', 'In_WaterBuffalo']].sum().reset_index()
    heat_data.columns = ["Species", "Gene Count"]
    plt.figure(figsize=(8, 5))
    sns.heatmap(heat_data[['Gene Count']].T, annot=True, fmt='d', cmap='YlGnBu', cbar=False,
                xticklabels=heat_data['Species'].tolist())
    plt.title("Total Gene Products Per Species")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gene_product_species_heatmap.png"))
    plt.close()

    # Save summary CSV
    summary_df.to_csv(os.path.join(output_dir, "gene_product_summary.csv"), index=False)
    print(f"[✓] Analysis complete! Files saved in: {output_dir}")


# === Paths to CSV files ===
base_path = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\7_Gene_extraction"

wildyak_path = f"{base_path}\\wildyak_CDS.csv"
takin_path = f"{base_path}\\takin_CDS.csv"
waterbuffalo_path = f"{base_path}\\waterbuffalo_CDS.csv"

# === Output folder for plots and CSVs ===
output_path = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\8_gene_visualization\gene_product_analysis_output"

# === Run the analysis ===
analyze_gene_product_overlap(
    wildyak_csv=wildyak_path,
    takin_csv=takin_path,
    buffalo_csv=waterbuffalo_path,
    output_dir=output_path
)
