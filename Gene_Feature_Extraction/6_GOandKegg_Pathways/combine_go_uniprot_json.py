import os
import json
import pandas as pd

# ---------------------------
# Update these paths as needed:
csv_path = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\6_GOandKegg_Pathways\compareCluster_GO\GO_compareCluster.csv"
uniprot_json_path = r"uniprot_results.json"  # Make sure this file is in the working directory or provide full path.
go_json_path = r"go_results.json"            # Make sure this file is in the working directory or provide full path.
output_csv_path = r"merged_enriched_GO.csv"
# ---------------------------

# Load the CSV file into a DataFrame.
df = pd.read_csv(csv_path)

# Load the UniProt JSON file.
with open(uniprot_json_path, "r") as f:
    uniprot_data = json.load(f)

# Load the GO JSON file.
with open(go_json_path, "r") as f:
    go_data = json.load(f)

def extract_go_info(go_id):
    """
    Look up a GO term's details in go_data.
    Returns a tuple: (name, definition, synonyms) if found,
    Otherwise ("N/A", "N/A", "").
    """
    if go_id in go_data:
        results = go_data[go_id].get("results", [])
        if results:
            entry = results[0]
            name = entry.get("name", "N/A")
            definition = entry.get("definition", {}).get("text", "N/A")
            syns = entry.get("synonyms", [])
            # Combine synonyms names separated by semicolon.
            synonyms = "; ".join([syn.get("name", "") for syn in syns]) if syns else ""
            return name, definition, synonyms
    return "N/A", "N/A", ""

def summarize_uniprot_info(gene_list):
    """
    Given a list of gene (protein) accessions,
    return a string summary containing the protein recommended name for each accession.
    If an entry contains an "error", note it.
    """
    summaries = []
    for gene in gene_list:
        gene = gene.strip()
        entry = uniprot_data.get(gene, {})
        if "error" in entry:
            summaries.append(f"{gene}: ERROR")
        else:
            # Extract recommended protein name if available.
            try:
                protein_name = (
                    entry.get("proteinDescription", {})
                        .get("recommendedName", {})
                        .get("fullName", {})
                        .get("value", "Unnamed")
                )
            except Exception:
                protein_name = "Unnamed"
            summaries.append(f"{gene}: {protein_name}")
    return " | ".join(summaries)

# Lists to store new columns.
go_names = []
go_definitions = []
go_synonyms = []
uniprot_summaries = []

# Process each row in the CSV.
for idx, row in df.iterrows():
    # Extract and clean the GO identifier from the "ID" column.
    raw_go_id = str(row["ID"])
    clean_go_id = raw_go_id.split("(")[0].strip()
    
    # Retrieve GO details.
    name, definition, synonyms = extract_go_info(clean_go_id)
    go_names.append(name)
    go_definitions.append(definition)
    go_synonyms.append(synonyms)
    
    # Process gene IDs from the "geneID" column (assumed to be "/" separated).
    raw_gene_ids = str(row["geneID"])
    gene_list = raw_gene_ids.split("/")
    up_summary = summarize_uniprot_info(gene_list)
    uniprot_summaries.append(up_summary)

# Add the new annotation columns to the DataFrame.
df["GO_Name"] = go_names
df["GO_Definition"] = go_definitions
df["GO_Synonyms"] = go_synonyms
df["UniProt_Summary"] = uniprot_summaries

# Write the enriched DataFrame to a new CSV file.
df.to_csv(output_csv_path, index=False)
print(f"Enriched CSV file saved as: {output_csv_path}")
