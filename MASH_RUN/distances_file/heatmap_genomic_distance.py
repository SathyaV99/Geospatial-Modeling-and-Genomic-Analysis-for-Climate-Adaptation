import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Mash distances file
df = pd.read_csv("highres_distances.tsv", sep="\t", header=None,
                 names=["Genome1", "Genome2", "Distance", "P-value", "Shared_kmers"])

# species name map
name_map = {
    "GCF_023091745.1_Takin1.1_genomic.fna": "Takin",
    "GCF_019923935.1_NDDB_SH_1_genomic.fna": "Water Buffalo",
    "GCF_027580195.1_NWIPB_WYAK_1.1_genomic.fna": "Wild Yak"
}

df["Genome1"] = df["Genome1"].map(name_map)
df["Genome2"] = df["Genome2"].map(name_map)

# present in the final matrix
genomes = ["Takin", "Water Buffalo", "Wild Yak"]
distance_matrix = pd.DataFrame(index=genomes, columns=genomes, data=0.0)

for _, row in df.iterrows():
    distance_matrix.at[row["Genome1"], row["Genome2"]] = row["Distance"]
    distance_matrix.at[row["Genome2"], row["Genome1"]] = row["Distance"]

for genome in genomes:
    distance_matrix.at[genome, genome] = 0.0

# Plot heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(distance_matrix, annot=True, fmt=".4f", cmap="viridis", square=True,
            linewidths=0.5, linecolor='gray', cbar_kws={'label': 'Mash Distance'})
plt.title("Genomic Distance Heatmap (Mash)")
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()
