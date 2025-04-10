import os
import pandas as pd
import requests
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

# Function to fetch protein information from UniProt
def fetch_uniprot_data(accession, timeout=5):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return accession, response.json()
    except Exception as e:
        return accession, {"error": str(e)}

# Function to fetch GO term details from QuickGO
def fetch_go_term(go_id, timeout=5):
    url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}"
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return go_id, response.json()
    except Exception as e:
        return go_id, {"error": str(e)}

if __name__ == '__main__':
    # Path to your CSV file:
    csv_path = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB\Gene_Feature_Extraction\6_GOandKegg_Pathways\compareCluster_GO\GO_compareCluster.csv"
    
    # Read the CSV file into a DataFrame.
    # Adjust the separator (sep) if necessary; here we try the default comma.
    df = pd.read_csv(csv_path)
    
    # For debugging, you can print the first few rows:
    # print(df.head())

    # Extract unique gene IDs from the 'geneID' column.
    gene_ids = set()
    for ids in df['geneID'].dropna():
        for gene in ids.split('/'):
            gene_ids.add(gene.strip())
    gene_ids = list(gene_ids)
    print(f"Found {len(gene_ids)} unique gene IDs.")

    # Extract unique GO IDs from the 'ID' column.
    # The ID column is in the format "GO:0004984(PANTHER)"; we remove the "(PANTHER)" part.
    go_ids = set()
    for go in df['ID'].dropna():
        go_id = go.split('(')[0].strip()  # "GO:0004984"
        go_ids.add(go_id)
    go_ids = list(go_ids)
    print(f"Found {len(go_ids)} unique GO term IDs.")

    # Dictionaries to store fetched API results
    uniprot_results = {}
    go_results = {}

    # Use ThreadPoolExecutor to fetch API data concurrently
    max_workers = 20  # Adjust the number of threads as needed

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit UniProt API calls concurrently for gene IDs
        future_to_gene = {executor.submit(fetch_uniprot_data, gene_id): gene_id for gene_id in gene_ids}
        for future in as_completed(future_to_gene):
            gene_id, data = future.result()
            uniprot_results[gene_id] = data

        # Submit QuickGO API calls concurrently for GO term IDs
        future_to_go = {executor.submit(fetch_go_term, go_id): go_id for go_id in go_ids}
        for future in as_completed(future_to_go):
            go_id, data = future.result()
            go_results[go_id] = data

    # Optionally, save the results to JSON files for later examination
    with open("uniprot_results.json", "w") as f:
        json.dump(uniprot_results, f, indent=2)
    with open("go_results.json", "w") as f:
        json.dump(go_results, f, indent=2)

    print("Completed fetching API data.")
