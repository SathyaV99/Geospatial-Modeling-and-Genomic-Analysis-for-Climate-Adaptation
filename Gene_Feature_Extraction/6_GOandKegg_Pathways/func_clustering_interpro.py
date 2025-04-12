import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import scipy.cluster.hierarchy as sch

def functional_clustering_interpro(csv_path, output_dir, n_clusters=4):
    """
    Performs functional clustering of InterPro domains across species.

    Parameters:
        csv_path (str): Path to the InterPro domain count CSV.
        output_dir (str): Output directory to save plots and clustered CSV.
        n_clusters (int): Number of clusters for KMeans (default=4).
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load the data
    df = pd.read_csv(csv_path)

    # Fill missing values with 0
    df.fillna(0, inplace=True)

    # Features = counts for each species (skip first column)
    features = df.columns[1:]
    X = df[features].values

    # Standardize the data
    X_scaled = StandardScaler().fit_transform(X)

    # === Hierarchical Clustering (Dendrogram) ===
    plt.figure(figsize=(14, 8))
    dendrogram = sch.dendrogram(
        sch.linkage(X_scaled, method='ward'),
        labels=df['InterPro_Desc'].values,
        leaf_rotation=90
    )
    plt.title("Hierarchical Clustering of InterPro Domains")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "interpro_hierarchical_clustering.png"))
    plt.close()

    # === KMeans Clustering ===
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(X_scaled)
    df['Cluster'] = clusters

    # Save clustered data
    clustered_path = os.path.join(output_dir, "interpro_clustered.csv")
    df.to_csv(clustered_path, index=False)

    # === Heatmap of Cluster Averages ===
    cluster_means = df.groupby('Cluster')[features].mean()
    plt.figure(figsize=(10, 6))
    sns.heatmap(cluster_means, annot=True, cmap="coolwarm", fmt=".0f")
    plt.title("Average InterPro Domain Counts per Cluster")
    plt.ylabel("Cluster")
    plt.tight_layout()
    heatmap_path = os.path.join(output_dir, "cluster_mean_heatmap.png")
    plt.savefig(heatmap_path)
    plt.close()

    print(f"[✓] Clustering complete!")
    print(f"📄 Clustered CSV: {clustered_path}")
    print(f"🧬 Dendrogram saved to: {output_dir}/interpro_hierarchical_clustering.png")
    print(f"📊 Heatmap saved to: {output_dir}/cluster_mean_heatmap.png")

    return df[['InterPro_Desc', 'Cluster']]

#_---------------------------
if __name__ == "__main__":
    csv_path = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\6_GOandKegg_Pathways\interpro_output\interpro_domain_comparison.csv"
    output_dir = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\6_GOandKegg_Pathways\interpro_output\functional_clustering_output"

    functional_clustering_interpro(csv_path, output_dir)
