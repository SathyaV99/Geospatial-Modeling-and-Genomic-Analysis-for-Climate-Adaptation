from Bio import SeqIO
import pandas as pd
import os

def extract_features_from_gbff(gbff_path, output_csv=None):
    """
    Extract annotated features from a GBFF (GenBank) file.

    Parameters:
        gbff_path (str): Path to the .gbff file
        output_csv (str, optional): Output CSV path. If None, auto-generates name.

    Returns:
        pd.DataFrame: DataFrame containing all features
    """
    features = []

    for record in SeqIO.parse(gbff_path, "genbank"):
        for feature in record.features:
            data = {
                "Contig": record.id,
                "Feature_Type": feature.type,
                "Start": int(feature.location.start),
                "End": int(feature.location.end),
                "Strand": feature.location.strand,
                "Locus_Tag": feature.qualifiers.get("locus_tag", [""])[0],
                "Gene": feature.qualifiers.get("gene", [""])[0],
                "Product": feature.qualifiers.get("product", [""])[0],
                "Protein_ID": feature.qualifiers.get("protein_id", [""])[0],
                "Translation": feature.qualifiers.get("translation", [""])[0] if "translation" in feature.qualifiers else "",
                "Note": feature.qualifiers.get("note", [""])[0] if "note" in feature.qualifiers else ""
            }
            features.append(data)

    df = pd.DataFrame(features)

    # Auto-generate output CSV name if not provided
    if not output_csv:
        base = os.path.splitext(os.path.basename(gbff_path))[0]
        output_csv = f"{base}_features.csv"

    # Ensure .csv extension is included
    if not output_csv.endswith(".csv"):
        output_csv += ".csv"

    df.to_csv(output_csv, index=False)
    print(f"[✓] Extracted features saved to: {output_csv}")

    return df

# List of GBFF files to process
gbff_files = [
    r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Gene_Feature_Extraction\takin_genomic.gbff",
    r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Gene_Feature_Extraction\waterbuffalo_genomic.gbff",
    r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Gene_Feature_Extraction\wildyak_genomic.gbff"
]

# Loop over and extract features from each GBFF file
for gbff in gbff_files:
    extract_features_from_gbff(gbff)
