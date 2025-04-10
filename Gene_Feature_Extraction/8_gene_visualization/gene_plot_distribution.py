import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_gene_distribution_barplot(summary_dict, output_path):
    # Create DataFrame
    df = pd.DataFrame(list(summary_dict.items()), columns=["Category", "Count"])
    total = df["Count"].sum()
    df["Percent"] = df["Count"] / total * 100

    # Bar Plot
    plt.figure(figsize=(12, 6))
    sns.set(style="whitegrid")
    ax = sns.barplot(data=df, x="Category", y="Count", palette="mako")

    # Annotate with count and %
    for index, row in df.iterrows():
        ax.text(index, row["Count"] + total * 0.01, f'{row["Count"]:,}\n({row["Percent"]:.1f}%)',
                color='black', ha='center', va='bottom', fontsize=9)

    plt.title("Gene Product Distribution Across Species", fontsize=14)
    plt.ylabel("Number of Gene Products")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"[✓] Barplot with % saved to: {output_path}")

# === Example Usage ===
summary = {
    "Shared (All 3)": 17822,
    "Unique to Wild Yak": 3542,
    "Unique to Takin": 2498,
    "Unique to Water Buffalo": 22098,
    "Only in Wild Yak (not Buffalo)": 13622,
    "Only in Water Buffalo (not Yak)": 30076,
    "Only in Takin (not Buffalo)": 15453,
}

output_barplot = "gene_product_distribution_barplot_with_percentage.png"
plot_gene_distribution_barplot(summary, output_barplot)
