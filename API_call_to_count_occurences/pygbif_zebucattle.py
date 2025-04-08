import requests
import pandas as pd

# GBIF taxon key for Bos indicus (Zebu cattle)
taxon_key = 2441023
limit = 300  # Number of records per request
offset = 0  # Starting offset
all_records = []

print("Fetching data from GBIF API for Zebu cattle (Bos indicus)...")

# Loop to fetch all occurrences using pagination
while True:
    url = f"https://api.gbif.org/v1/occurrence/search?taxonKey={taxon_key}&limit={limit}&offset={offset}"
    response = requests.get(url)
    data = response.json()

    # Extract relevant fields
    records = [
        {
            "species": rec.get("scientificName", "Unknown"),
            "latitude": rec.get("decimalLatitude", "No Coordinates"),
            "longitude": rec.get("decimalLongitude", "No Coordinates"),
            "country": rec.get("country", "Unknown"),
            "year": rec.get("year", "Unknown"),
            "dataset": rec.get("datasetName", "Unknown"),
            "basisOfRecord": rec.get("basisOfRecord", "Unknown")
        }
        for rec in data.get("results", [])
    ]

    # Add to the main list
    all_records.extend(records)

    # Break the loop if there are no more records
    if len(data.get("results", [])) < limit:
        break  # Stop when no more pages exist

    # Increase offset for the next batch
    offset += limit
    print(f"Retrieved {len(all_records)} records so far...")

# Convert to DataFrame
df = pd.DataFrame(all_records)

# Save to CSV
csv_filename = "zebu_cattle_gbif_occurrences_full.csv"
df.to_csv(csv_filename, index=False)

print(f"Total {len(df)} records downloaded. Data saved to {csv_filename}.")
