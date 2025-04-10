import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def analyze_gene_family_variation(input_file, output_dir, species_order=["wild_yak", "takin", "water_buffalo"], cv_threshold=0.5):
    """
    Analyzes gene family counts (orthogroups) across species to detect expansion/contraction signals.
    
    The input_file should be a TSV file with columns:
         FamilyID, wild_yak, takin, water_buffalo
    (Order must match the species tree: wild_yak, takin, water_buffalo)
    
    The function computes for each gene family:
      - Mean, standard deviation, CV (std/mean), max, min, and range.
      - Generates a histogram of the CV values.
      - Generates a scatter plot of wild_yak vs. takin counts.
      - Filters gene families with a coefficient of variation (CV) above cv_threshold.
    
    Outputs:
      - A TSV file with computed statistics for every gene family.
      - A TSV file listing gene families with high variability.
      - A histogram plot of CV.
      - A scatter plot (wild_yak vs. takin).
      
    Parameters:
      input_file (str): Path to the input TSV file.
      output_dir (str): Directory to save output files.
      species_order (list): List of species column names. (Default: ["wild_yak", "takin", "water_buffalo"])
      cv_threshold (float): Minimum CV value to flag high variability (default 0.5).
    
    Returns:
      DataFrame: DataFrame with per-gene-family statistics.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Load the input file (TSV) with proper separator
    df = pd.read_csv(input_file, sep="\t")
    
    # Ensure species columns are numeric
    for sp in species_order:
        df[sp] = pd.to_numeric(df[sp], errors="coerce").fillna(0)
    
    # Compute statistics for each gene family
    df['mean'] = df[species_order].mean(axis=1)
    df['std'] = df[species_order].std(axis=1)
    # Avoid division by zero: if mean==0, set CV to 0; else CV = std/mean
    df['cv'] = np.where(df['mean'] > 0, df['std'] / df['mean'], 0)
    df['max'] = df[species_order].max(axis=1)
    df['min'] = df[species_order].min(axis=1)
    df['range'] = df['max'] - df['min']
    
    # Save the full statistics table
    stats_file = os.path.join(output_dir, "gene_family_stats.tsv")
    df.to_csv(stats_file, sep="\t", index=False)
    print(f"[✓] Gene family stats saved to: {stats_file}")
    
    # Plot histogram of coefficient of variation (CV)
    plt.figure(figsize=(8, 6))
    sns.histplot(df['cv'], bins=30, kde=True, color="steelblue")
    plt.xlabel("Coefficient of Variation (CV)")
    plt.title("Distribution of Gene Family CV Across Species")
    hist_cv_path = os.path.join(output_dir, "gene_family_cv_histogram.png")
    plt.tight_layout()
    plt.savefig(hist_cv_path)
    plt.close()
    print(f"[✓] CV Histogram saved to: {hist_cv_path}")
    
    # Plot scatter plot: wild_yak vs takin counts
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=df, x="wild_yak", y="takin", alpha=0.6, color="darkorange")
    plt.xlabel("Wild Yak Gene Count")
    plt.ylabel("Takin Gene Count")
    plt.title("Wild Yak vs. Takin Gene Family Counts")
    scatter_path = os.path.join(output_dir, "scatter_wild_yak_vs_takin.png")
    plt.tight_layout()
    plt.savefig(scatter_path)
    plt.close()
    print(f"[✓] Scatter plot saved to: {scatter_path}")
    
    # Identify gene families with high variability (CV above threshold)
    high_variation = df[df['cv'] >= cv_threshold]
    print(f"[✓] Found {len(high_variation)} gene families with CV >= {cv_threshold}")
    high_var_file = os.path.join(output_dir, "high_variation_gene_families.tsv")
    high_variation.to_csv(high_var_file, sep="\t", index=False)
    print(f"[✓] High variability gene families saved to: {high_var_file}")
    
    return df

# === Example Usage ===
input_file = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\cafe_input_gene_counts_final_fixed_final3.tsv"
output_dir = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\gene_family_variation_output1"

df_stats = analyze_gene_family_variation(input_file, output_dir, species_order=["wild_yak", "takin", "water_buffalo"], cv_threshold=0.5)
