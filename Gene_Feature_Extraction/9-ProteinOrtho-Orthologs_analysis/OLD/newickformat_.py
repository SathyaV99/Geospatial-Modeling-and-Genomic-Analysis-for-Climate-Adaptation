def save_newick_tree(file_path):
    newick_tree = "((wild_yak:0.02,takin:0.02):0.03,water_buffalo:0.05);"
    with open(file_path, 'w') as f:
        f.write(newick_tree)
    print(f"[✓] Newick tree saved to: {file_path}")

# === Example path to save ===
output_path = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\9-ProteinOrtho-Orthologs_analysis\cafe_species_tree.nwk"
save_newick_tree(output_path)
