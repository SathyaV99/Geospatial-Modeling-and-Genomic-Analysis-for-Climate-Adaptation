import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os

def load_and_prepare_proteinortho(tsv_path):

    # Read TSV with comment #
    df = pd.read_csv(tsv_path, sep="\t", header=None, comment="#")
    
    headers = ["FamilyID", "water_buffalo", "wild_yak", "takin"]
    df = df.iloc[:, :4]
    df.columns = headers
    return df

def create_presence_matrix(df):

    presence_df = df.copy()
    for col in ["water_buffalo", "wild_yak", "takin"]:
        # 1 or 0
        presence_df[col] = presence_df[col].notna() & (presence_df[col] != "")
        presence_df[col] = presence_df[col].astype(int)
    return presence_df

def summarize_orthogroups(presence_df):

    total = len(presence_df)
    core = len(presence_df[(presence_df["water_buffalo"]==1) & 
                            (presence_df["wild_yak"]==1) & 
                            (presence_df["takin"]==1)])
    unique_buffalo = len(presence_df[(presence_df["water_buffalo"]==1) & 
                                     (presence_df["wild_yak"]==0) & 
                                     (presence_df["takin"]==0)])
    unique_yak = len(presence_df[(presence_df["water_buffalo"]==0) & 
                                 (presence_df["wild_yak"]==1) & 
                                 (presence_df["takin"]==0)])
    unique_takin = len(presence_df[(presence_df["water_buffalo"]==0) & 
                                   (presence_df["wild_yak"]==0) & 
                                   (presence_df["takin"]==1)])
    
    # Shared by 2 species:
    shared_buffalo_yak = len(presence_df[(presence_df["water_buffalo"]==1) & 
                                         (presence_df["wild_yak"]==1) & 
                                         (presence_df["takin"]==0)])
    shared_buffalo_takin = len(presence_df[(presence_df["water_buffalo"]==1) & 
                                           (presence_df["wild_yak"]==0) & 
                                           (presence_df["takin"]==1)])
    shared_yak_takin = len(presence_df[(presence_df["water_buffalo"]==0) & 
                                       (presence_df["wild_yak"]==1) & 
                                       (presence_df["takin"]==1)])
    
    summary = {
        "Total orthogroups": total,
        "Core (all 3 present)": core,
        "Unique to Water Buffalo": unique_buffalo,
        "Unique to Wild Yak": unique_yak,
        "Unique to Takin": unique_takin,
        "Shared by WaterBuffalo & WildYak": shared_buffalo_yak,
        "Shared by WaterBuffalo & Takin": shared_buffalo_takin,
        "Shared by WildYak & Takin": shared_yak_takin,
    }
    
    return summary

def plot_venn(presence_df, output_path):
    """
    Plots a Venn diagram of orthogroup sharing among the three species.
    """
    buffalo_set = set(presence_df.loc[presence_df["water_buffalo"]==1, "FamilyID"])
    yak_set = set(presence_df.loc[presence_df["wild_yak"]==1, "FamilyID"])
    takin_set = set(presence_df.loc[presence_df["takin"]==1, "FamilyID"])
    
    plt.figure(figsize=(8,6))
    venn3([buffalo_set, yak_set, takin_set], set_labels=("Water Buffalo", "Wild Yak", "Takin"))
    plt.title("Orthogroup Sharing Among Species")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f":D Venn diagram saved to: {output_path}")

def plot_bar_summary(summary_dict, output_path):
    """
    Plots a barplot with both absolute counts and percentages of gene family sharing.
    """
    df_summary = pd.DataFrame(list(summary_dict.items()), columns=["Category", "Count"])
    total = df_summary["Count"].sum()
    df_summary["Percent"] = df_summary["Count"] / total * 100
    
    plt.figure(figsize=(10,6))
    ax = sns.barplot(data=df_summary, x="Category", y="Count", palette="viridis")
    plt.xticks(rotation=45, ha="right")
    plt.title("Orthogroup Distribution Summary")
    
    # Annotate bars with count and percent
    for idx, row in df_summary.iterrows():
        ax.text(idx, row["Count"]+total*0.005, f'{row["Count"]}\n({row["Percent"]:.1f}%)',
                ha="center", va="bottom", fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f":D Barplot summary saved to: {output_path}")

def plot_pca(presence_df, output_path):
    """
    Performs PCA on the presence/absence matrix and plots a 2D scatter plot.
    Each row is an orthogroup; though this matrix is binary, PCA helps visualize variance.
    """
    X = presence_df[["water_buffalo", "wild_yak", "takin"]].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(X_scaled)
    
    pca_df = pd.DataFrame(data=principalComponents, columns=["PC1", "PC2"])
    pca_df["FamilyID"] = presence_df["FamilyID"]
    
    plt.figure(figsize=(8,6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", s=10, color="darkblue")
    plt.title("PCA of Orthogroup Presence/Absence")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f":D PCA plot saved to: {output_path}")

if __name__ == "__main__":
    # FILE PATHS
    proteinortho_tsv = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\proteinortho\myproject.proteinortho.tsv"
    output_folder = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\orthogroup_visualization"
    os.makedirs(output_folder, exist_ok=True)

    df_ortho = load_and_prepare_proteinortho(proteinortho_tsv)
    
    # Create a presence/absence matrix
    presence_df = create_presence_matrix(df_ortho)
    
    # Generate summary statistics
    summary = summarize_orthogroups(presence_df)
    print("Summary Statistics:")
    for key, value in summary.items():
        print(f"{key}: {value}")
    
    # SAVE MATRIX
    presence_matrix_path = os.path.join(output_folder, "orthogroup_presence_matrix.csv")
    presence_df.to_csv(presence_matrix_path, index=False)
    print(f"[✓] Presence matrix saved to: {presence_matrix_path}")
    
    # VENN DIAgram
    venn_output = os.path.join(output_folder, "venn_orthogroups.png")
    plot_venn(presence_df, venn_output)
    
    # summary
    bar_output = os.path.join(output_folder, "barplot_orthogroup_summary.png")
    plot_bar_summary(summary, bar_output)
    
    # pca su
    pca_output = os.path.join(output_folder, "pca_orthogroup_presence.png")
    plot_pca(presence_df, pca_output)
