import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns
import os

def load_presence_matrix(matrix_csv):
    """
    Load the presence matrix CSV.
    If the file is read in as a single column (with literal "\t"),
    this function will split it into separate columns.
    """
    # Attempt to read using tab as delimiter
    try:
        df = pd.read_csv(matrix_csv, sep="\t")
    except Exception as e:
        print("Error reading file with sep='\\t':", e)
        df = pd.read_csv(matrix_csv)
    
    # If only one column is read, check if it contains literal "\t" and split it.
    if df.shape[1] == 1:
        # Split the first column on the literal "\t"
        df = df[df.columns[0]].str.split("\\t", expand=True)
    
    # Return the dataframe
    return df

def visualize_presence_matrix_pca_and_bar(matrix_csv, output_pca, output_bar, adaptive_keywords=None):
    """
    Reads a gene family presence/absence matrix CSV file,
    runs PCA on the presence columns, and generates:
      1. A scatter plot of the first two principal components,
         highlighting adaptive gene families if specified.
      2. A barplot showing both absolute counts and percentages per species.
    
    The CSV file should contain an identifier column (e.g., "FamilyID" or "Gene_Product")
    and numeric columns for the species named 'wild_yak', 'takin', 'water_buffalo' (case-insensitive).
    
    Parameters:
        matrix_csv (str): Path to the CSV file.
        output_pca (str): Path (including filename) to save the PCA scatter plot.
        output_bar (str): Path (including filename) to save the barplot.
        adaptive_keywords (list of str): List of keywords to flag adaptive genes.
    """
    # Load the CSV file robustly
    df = load_presence_matrix(matrix_csv)
    
    # Normalize column names: strip extra spaces and convert to lower case.
    df.columns = df.columns.astype(str).str.strip().str.lower()
    
    # Determine the identifier column: check for 'gene_product' or 'familyid', else use first column.
    if 'gene_product' in df.columns:
        id_col = 'gene_product'
    elif 'familyid' in df.columns:
        id_col = 'familyid'
    else:
        id_col = df.columns[0]
    
    # Expected species columns (lowercase)
    expected_cols = ['wild_yak', 'takin', 'water_buffalo']
    for col in expected_cols:
        if col not in df.columns:
            raise KeyError(f"Expected column '{col}' not found. Available columns: {df.columns.tolist()}")
    
    # Ensure the presence columns are numeric (0 or 1)
    for col in expected_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
    
    # ---------- PCA Visualization ----------
    X = df[expected_cols].values
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)
    df['pc1'] = pcs[:, 0]
    df['pc2'] = pcs[:, 1]
    
    plt.figure(figsize=(10, 8))
    plt.scatter(df['pc1'], df['pc2'], color='gray', alpha=0.5, label='Gene Families')
    
    # Highlight adaptive gene families if adaptive_keywords is provided
    if adaptive_keywords:
        mask = df[id_col].astype(str).str.lower().apply(
            lambda x: any(keyword.lower() in x for keyword in adaptive_keywords)
        )
        plt.scatter(df.loc[mask, 'pc1'], df.loc[mask, 'pc2'], color='red', alpha=0.8, label='Adaptive Genes')
    
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA of Gene Family Presence/Absence')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_pca)
    plt.close()
    print(f"[✓] PCA plot saved to {output_pca}")
    
# ---------- Barplot Visualization (Taller) ----------
    counts = df[expected_cols].sum()
    total = counts.sum()
    percents = counts / total * 100
    bar_df = pd.DataFrame({
        'Species': counts.index,
        'Count': counts.values,
        'Percentage': percents.values
    })

# Increase the height by adjusting figsize=(width, height)
    plt.figure(figsize=(8, 15))  # 6 wide, 10 tall
    ax = sns.barplot(data=bar_df, x='Species', y='Count', palette="mako")

    plt.title('Total Gene Family Counts Across Species')
    plt.ylabel('Absolute Count')
    plt.xlabel('')

# Annotate each bar with count and percentage
    for idx, row in bar_df.iterrows():
        ax.text(
            idx,
            row['Count'] + total * 0.01,
            f"{int(row['Count']):,}\n({row['Percentage']:.1f}%)",
            color='black',
            ha='center',
            va='bottom',
            fontsize=10
        )

    plt.tight_layout()
    plt.savefig(output_bar)  # your existing path
    plt.close()
    print(f"[✓] Barplot saved to {output_bar}")

    
    return df

# === Example usage ===
# Update these paths to match your file locations
matrix_csv = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\cafe_input_gene_counts_final_fixed_final3.tsv"
output_pca = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\gene_presence_pca.png"
output_bar = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\gene_presence_barplot.png"
adaptive_keywords = ["hemoglobin", "cold", "metab", "hypoxia", "oxygen", "heat", "stress", "HIF"]

df_result = visualize_presence_matrix_pca_and_bar(matrix_csv, output_pca, output_bar, adaptive_keywords)
