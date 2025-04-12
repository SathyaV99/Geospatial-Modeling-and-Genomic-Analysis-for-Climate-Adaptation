import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns
import os

def load_presence_matrix(matrix_csv):
    try:
        df = pd.read_csv(matrix_csv, sep="\t")
    except Exception as e:
        print("Error reading file with sep='\\t':", e)
        df = pd.read_csv(matrix_csv)
 
    if df.shape[1] == 1:
        # Split the first column on the literal "\t"
        df = df[df.columns[0]].str.split("\\t", expand=True)
    return df

def visualize_presence_matrix_pca_and_bar(matrix_csv, output_pca, output_bar, adaptive_keywords=None):
    df = load_presence_matrix(matrix_csv)
    
    # Normalize column names
    df.columns = df.columns.astype(str).str.strip().str.lower()
    
    if 'gene_product' in df.columns:
        id_col = 'gene_product'
    elif 'familyid' in df.columns:
        id_col = 'familyid'
    else:
        id_col = df.columns[0]
    
    # Expected species columns
    expected_cols = ['wild_yak', 'takin', 'water_buffalo']
    for col in expected_cols:
        if col not in df.columns:
            raise KeyError(f"Expected column '{col}' not found. Available columns: {df.columns.tolist()}")
    
    #  1 or 0 for spec name
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
    print(f":D PCA plot saved to {output_pca}")
    
# ---------- Barplot Visualization (Taller) ----------
    counts = df[expected_cols].sum()
    total = counts.sum()
    percents = counts / total * 100
    bar_df = pd.DataFrame({
        'Species': counts.index,
        'Count': counts.values,
        'Percentage': percents.values
    })

    plt.figure(figsize=(8, 15))  # 6 wide, 10 tall
    ax = sns.barplot(data=bar_df, x='Species', y='Count', palette="mako")

    plt.title('Total Gene Family Counts Across Species')
    plt.ylabel('Absolute Count')
    plt.xlabel('')

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
    plt.savefig(output_bar) 
    plt.close()
    print(f":D Barplot saved to {output_bar}")

    
    return df

# PATHS
matrix_csv = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\cafe_input_gene_counts_final_fixed_final3.tsv"
output_pca = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\gene_presence_pca.png"
output_bar = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\gene_presence_barplot.png"
adaptive_keywords = ["hemoglobin", "cold", "metab", "hypoxia", "oxygen", "heat", "stress", "HIF"]

df_result = visualize_presence_matrix_pca_and_bar(matrix_csv, output_pca, output_bar, adaptive_keywords)
