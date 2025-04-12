import pandas as pd
import matplotlib.pyplot as plt
import os

plt.style.use('grayscale')

def creative_bar_grid(csv_path, species_name):
    # Load data from CSV
    df = pd.read_csv(csv_path)

    # 1. Feature Type Counts (top 10)
    feature_counts = df['Feature_Type'].value_counts().head(10)
    
    # 2. Top 10 Genes
    top_genes = df['Gene'].value_counts().head(10)
    
    # 3. Top 10 Products (ignoring missing product fields)
    if 'Product' in df.columns:
        top_products = df['Product'].dropna().value_counts().head(10)
    else:
        top_products = pd.Series(dtype=float)
    
    # 4. Strand Distribution (top 10)
    strand_counts = df['Strand'].value_counts().head(10)
    
    # truncate labels longer than 25 characters
    def truncate_labels(labels):
        return [label[:25] if len(label) > 25 else label for label in labels]
    
    # Create a figure with a 2x2 grid of subplots
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 12))
    fig.suptitle(f"{species_name} - Genomic Features Bar Plots", fontsize=16)
    
    # ---------------------------------------------------------
    # Plot 1: Feature Type Counts
    # ---------------------------------------------------------
    ax = axes[0][0]
    x_vals = list(range(len(feature_counts)))
    truncated_feature_labels = truncate_labels(feature_counts.index.astype(str))
    bars = ax.bar(x_vals, feature_counts.values, edgecolor='black', color='white')
    for bar in bars:
        bar.set_hatch('//')
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f'{int(height)}',
            ha='center',
            va='bottom',
            fontsize=8
        )
    ax.set_title("Feature Type Counts")
    ax.set_xlabel("Feature Type")
    ax.set_ylabel("Count")
    ax.set_xticks(x_vals)
    ax.set_xticklabels(truncated_feature_labels, rotation=45)
    
    # ---------------------------------------------------------
    # Plot 2: Top 10 Genes
    # ---------------------------------------------------------
    ax = axes[0][1]
    x_vals = list(range(len(top_genes)))
    truncated_gene_labels = truncate_labels(top_genes.index.astype(str))
    bars = ax.bar(x_vals, top_genes.values, edgecolor='black', color='white')
    for bar in bars:
        bar.set_hatch('xx')
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f'{int(height)}',
            ha='center',
            va='bottom',
            fontsize=8
        )
    ax.set_title("Top 10 Genes")
    ax.set_xlabel("Gene")
    ax.set_ylabel("Frequency")
    ax.set_xticks(x_vals)
    ax.set_xticklabels(truncated_gene_labels, rotation=45)
    
    # ---------------------------------------------------------
    # Plot 3: Top 10 Products
    # ---------------------------------------------------------
    ax = axes[1][0]
    if not top_products.empty:
        x_vals = list(range(len(top_products)))
        truncated_product_labels = truncate_labels(top_products.index.astype(str))
        bars = ax.bar(x_vals, top_products.values, edgecolor='black', color='white')
        for bar in bars:
            bar.set_hatch('||')
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height,
                f'{int(height)}',
                ha='center',
                va='bottom',
                fontsize=8
            )
        ax.set_title("Top 10 Products")
        ax.set_xlabel("Product")
        ax.set_ylabel("Frequency")
        ax.set_xticks(x_vals)
        ax.set_xticklabels(truncated_product_labels, rotation=45)
    else:
        ax.text(0.5, 0.5, "No Product data available",
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes)
        ax.set_axis_off()
    
    # ---------------------------------------------------------
    # Plot 4: Strand Distribution
    # ---------------------------------------------------------
    ax = axes[1][1]
    x_vals = list(range(len(strand_counts)))
    truncated_strand_labels = truncate_labels(strand_counts.index.astype(str))
    bars = ax.bar(x_vals, strand_counts.values, edgecolor='black', color='white')
    for bar in bars:
        bar.set_hatch('\\\\')
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f'{int(height)}',
            ha='center',
            va='bottom',
            fontsize=8
        )
    ax.set_title("Strand Distribution")
    ax.set_xlabel("Strand")
    ax.set_ylabel("Count")
    ax.set_xticks(x_vals)
    ax.set_xticklabels(truncated_strand_labels, rotation=45)
    
    # Adjust layout and save figure
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    output_path = os.path.splitext(csv_path)[0] + f"_{species_name}_bargrid.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

# FILE PATH
files = {
    "Takin": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\1_genomic_feature_extraction\takin_genomic_features.csv",
    "WaterBuffalo": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\1_genomic_feature_extraction\waterbuffalo_genomic_features.csv",
    "WildYak": r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\1_genomic_feature_extraction\wildyak_genomic_features.csv"
}

# save the bar grid image for each species
for species, path in files.items():
    output_file = creative_bar_grid(path, species)
    print(f"Bar grid for {species} saved to: {output_file}")
