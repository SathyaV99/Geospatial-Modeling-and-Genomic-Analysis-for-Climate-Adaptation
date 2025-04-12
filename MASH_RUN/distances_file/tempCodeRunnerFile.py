import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load THE Mash results
df = pd.read_csv("highres_distances.tsv", sep="\t", header=None,
                 names=["Genome1", "Genome2", "Distance", "p", "Shared"])

# Pivot into matrix
pivot = df.pivot(index="Genome1", columns="Genome2", values="Distance")
pivot = pivot.fillna(pivot.T)  # Fill symmetric values

# Plot heatmap
sns.heatmap(pivot, cmap="viridis", annot=True)
plt.title("Genomic Distance Heatmap (Mash)")
plt.show()
