import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import seaborn as sns
import os

def visualize_orthogroup_distributions(proteinortho_path, output_dir):
    # Setup
    os.makedirs(output_dir, exist_ok=True)
    headers = ["Orthogroup_ID", "Connectivity", "Species_Count", "Water_Buffalo", "Wild_Yak", "Takin"]
    df = pd.read_csv(proteinortho_path, sep="\t", header=None, names=headers, comment="#")

    # Convert to presence/absence matrix
    df['In_Buffalo'] = df['Water_Buffalo'].notna()
    df['In_Yak'] = df['Wild_Yak'].notna()
    df['In_Takin'] = df['Takin'].notna()

    # Save full presence/absence matrix
    presence_df = df[['Orthogroup_ID', 'In_Buffalo', 'In_Yak', 'In_Takin']]
    presence_df.to_csv(os.path.join(output_dir, "orthogroup_presence_matrix.csv"), index=False)

    # === 1. VENN DIAGRAM ===
    buffalo_set = set(df[df['In_Buffalo']]['Orthogroup_ID'])
    yak_set = set(df[df['In_Yak']]['Orthogroup_ID'])
    takin_set = set(df[df['In_Takin']]['Orthogroup_ID'])

    plt.figure(figsize=(8, 6))
    venn3([buffalo_set, yak_set, takin_set], set_labels=('Water Buffalo', 'Wild Yak', 'Takin'))
    plt.title("Orthogroup Venn Diagram")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "venn_orthogroups.png"))
    plt.close()

    # === 2. CORE vs ACCESSORY BARPLOT ===
    total = len(df)
    core = len(df[df['In_Buffalo'] & df['In_Yak'] & df['In_Takin']])
    buffalo_only = len(buffalo_set - (yak_set | takin_set))
    yak_only = len(yak_set - (buffalo_set | takin_set))
    takin_only = len(takin_set - (buffalo_set | yak_set))
    shared_some = total - core - buffalo_only - yak_only - takin_only

    bar_data = pd.DataFrame({
        'Category': ['Core', 'Shared (not all)', 'Unique to Buffalo', 'Unique to Yak', 'Unique to Takin'],
        'Count': [core, shared_some, buffalo_only, yak_only, takin_only]
    })

    plt.figure(figsize=(10, 6))
    sns.barplot(data=bar_data, x="Category", y="Count", palette="viridis")
    plt.title("Orthogroup Distribution Across Species")
    plt.ylabel("Number of Orthogroups")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "barplot_core_vs_accessory.png"))
    plt.close()

    # === 3. HEATMAP ===
    heat_df = presence_df.set_index("Orthogroup_ID").astype(int)
    plt.figure(figsize=(12, 6))
    sns.heatmap(heat_df.head(100), cmap="coolwarm", cbar=True)  # Show top 100 for clarity
    plt.title("Orthogroup Presence/Absence (First 100)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "orthogroup_heatmap.png"))
    plt.close()

    print(":DVenn diagram, barplot, and heatmap saved to:", output_dir)

proteinortho_file = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\proteinortho\myproject.proteinortho.tsv"
output_folder = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\orthogroup_visualization"

visualize_orthogroup_distributions(proteinortho_file, output_folder)
