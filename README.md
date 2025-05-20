# Genomic and Geographic Adaptations to Climate Change: A Comparative Study of Wild Yak, Takin, and High-Altitude Bovids 

## Wild Yak faces substantial habitat loss and altitudinal displacement by 2050, driven by climate change and limited adaptive genomic traits.

### Computationally intensive analyses like InterProScan and BLAST averaged 2–4 days each for a single species to complete!

This project investigates how high-altitude bovids-Wild Yak, Takin, and Water Buffalo-adapt to climate change. It combines species distribution modeling and comparative genomics to find out:

* Where their habitats are now and where they'll likely shift in the future

* What genes and protein functions help them survive in extreme environments

It uses Python and R with machine learning (Random Forest), spatial analysis (centroid tracking, jittering), and genomic tools like Mash, InterProScan, ProteinOrtho, and GO enrichment.

![image](https://github.com/user-attachments/assets/dfb5b206-fa64-4ac6-9397-ed88ebc2df1b)

**Key finding for SDM:**


* Wild yak habitats are shrinking and moving uphill. Takin habitats are expanding. These results help guide future conservation efforts.

![image](https://github.com/user-attachments/assets/668e54af-9175-413b-8cba-680e2d794254)

**Key finding for Genomic Analysis:**

* Wild yaks have fewer heat shock genes, and due to their thick fur, they struggle to regulate body heat. They are adapted to cold climates at elevations around 3,000 feet and cannot tolerate warmer temperatures.

<img width="1012" alt="immune_genes" src="https://github.com/user-attachments/assets/64fc52cd-2be0-4b08-83bd-89689bc10c3c" />

# 🌍 Species Distribution Modeling (SDM)

This part of the project models current and future habitat suitability for **Wild Yak** and **Takin** using geospatial and climate data. It applies machine learning to predict where these animals can survive based on environmental conditions.

---

### 🔸 Objective

Predict species range shifts from **2009 to 2050** using environmental variables and occurrence records.

---

## 📦 Data Sources

## 🌦️ Climate and Environmental Data

High-resolution environmental variables used to model habitat suitability for Wild Yak and Takin.

### Data Sources

- **TerraClimate (2009–2024)** – [https://www.climatologylab.org/terraclimate.html](https://www.climatologylab.org/terraclimate.html)
- **WorldClim 2050 SSP245 & SSP585** – [https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html)
- **Google Earth Engine DEM** – [https://developers.google.com/earth-engine/datasets](https://developers.google.com/earth-engine/datasets)
- **Natural Earth Landmask** – [https://www.naturalearthdata.com](https://www.naturalearthdata.com)

---

### Environmental Data Summary

- **Precipitation**: TerraClimate `.nc` files (2009–2024), annual sum, stacked.
- **Min Temperature**: TerraClimate `.nc` files (2009–2024), annual mean.
- **Max Temperature**: TerraClimate `.nc` files (2009–2024), annual mean.
- **Future Climate**: WorldClim SSP245 and SSP585 for 2050.
- **Elevation**: Merged `.tif` from Earth Engine, resampled with GDAL warp.
- **Landmask**: Rasterized from Natural Earth shapefile.

---

### Processing Steps

- Download monthly ppt, tmin, tmax NetCDF files.
- Aggregate annual values using `xarray`.
- Merge and resample elevation `.tif` files using GDAL warp.
- Rasterize Asia land shapefile to create landmask.
- Align all layers to the same spatial grid.

## 📍 Occurrence Data

Species presence data used for SDM modeling.

- **Species**:
  - Wild Yak: 366 records
  - Takin: 692 records

- **Columns**:
  - Longitude, Latitude
  - Station Name, Climate ID, Date/Time, Year, Month, Day
  - Max/Min/Mean Temp (°C)
  - Heat/Cool Degree Days (°C)
  - Total Rain (mm), Total Snow (cm), Total Precip (mm)
  - Snow on Ground (cm)
  - Wind Gust Direction and Speed
  - Data Quality Flags

- Data cleaned and spatially jittered.
- Combined with environmental layers for model input.

---

### 🛠️ Methodology

#### 1. **Data Preprocessing**
- Downloaded and cleaned species presence data (lat/lon, date).
- Applied **spatial jittering** to reduce location bias:
  - Wild Yak: 10 synthetic points per record
  - Takin: 2 synthetic points per record
- Climate variables:
  - **Total Precipitation**
  - **Minimum Temperature**
  - **Maximum Temperature**
  - **Elevation** (resampled to climate resolution)
- Pseudo-absence points generated randomly.

#### 2. **Modeling**
- **Algorithm Used**: Random Forest Classifier (`scikit-learn`)
- **Training/Test Split**: 70/30
- **Evaluation Metrics**: ROC-AUC, confusion matrix
- **Best ROC-AUC**:
  - Wild Yak: **0.999**
  - Takin: **0.98+**

#### 3. **Prediction & Mapping**
- Suitability scores from **0 to 1** generated for each year (2009–2024).
- Future projections mapped using SSP245 and SSP585 climate scenarios (2050).
- Threshold (0.5) used to classify presence/absence.
- Habitat centroids calculated annually to track spatial shifts.

---

### 📈 Key Results

| Species     | Trend                           | Elevation Shift     | Centroid Movement |
|-------------|----------------------------------|----------------------|--------------------|
| Wild Yak    | Habitat **shrinks** by 2050      | 4750m → ~4810m       | NW by ~110 km      |
| Takin       | Habitat **expands** by 2050      | Increase expected    | W by ~121 km        |

---

### 🗂️ File Structure

```
SDM/
├── Code/
│   └── SDM_Final.ipynb        # Core modeling notebook
├── Data/
│   ├── takin_Final_cleaned.xls
│   ├── wild_yak_Final_cleaned.xls
│   └── elevation_resampled_to_climate.tif
├── Output/
│   ├── sdm_takin/
│   │   ├── suitability_map_20XX.png, .npy
│   │   ├── centroid_shifts_takin.csv
│   │   └── takin_suitability_area_trend.png
│   └── sdm_yak/
│       ├── suitability_map_20XX.png, .npy
│       ├── centroid_shifts.csv
│       └── yak_suitability_area_trend_final.png
```

---
### ▶️ How to Run

```bash
cd SDM/Code
jupyter notebook SDM_Final.ipynb
```

Make sure to have the following Python packages installed:
```bash
pip install scikit-learn rasterio xarray numpy pandas matplotlib
```

---

# 🐝 Genomic Analysis

This part of the project investigates **genetic adaptations** of Wild Yak, Takin, and Water Buffalo by comparing their full genomes, protein domains, and gene families.

---

### 🔸 Objective

Identify genetic traits linked to high-altitude survival using comparative genomics and functional annotations.

---

### 📦 Data Sources

To perform a comprehensive comparative genomics analysis of high-altitude bovids, we curated and processed full genome assemblies and annotations for three target species: **Wild Yak** (*Bos mutus*), **Takin** (*Budorcas taxicolor*), and **Water Buffalo** (*Bubalus bubalis*). These datasets were downloaded from the NCBI Assembly and GenBank/RefSeq repositories, using the most recent high-quality assemblies.

### 🔗 Genome Assembly Links

- [Water Buffalo](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=89462)
- [Indicine Cattle](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9915)
- [Wild Yak](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=72004)
- [Takin](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=37181)

---

### 📁 3.2.1 Genome File Types and Attributes

#### i) Original Datasets

Each genome dataset from NCBI contains multiple standard annotation files:

| File Type              | Format | Description                                 |
|------------------------|--------|---------------------------------------------|
| `.genomic.fna`         | FASTA  | Whole-genome nucleotide sequence            |
| `cds_from_genomic.fna` | FASTA  | Coding sequences (CDS)                      |
| `genomic.gbff`         | GBFF   | GenBank flat file with annotations          |
| `genomic.gff / .gtf`   | GFF/GTF| Gene feature coordinates                    |
| `*.faa` (derived)      | FASTA  | Translated protein sequences                |

**Genome file sizes (compressed/uncompressed):**

| Species         | Compressed (GB) | Uncompressed (GB) |
|-----------------|------------------|--------------------|
| Wild Yak        | 3.60 GB          | 12.10 GB           |
| Takin           | 3.73 GB          | 12.60 GB           |
| Water Buffalo   | 3.63 GB          | 12.70 GB           |

---

#### ii) Derived Datasets

To support analysis, original files were parsed and converted into structured, readable formats.

##### a) Genomic Features Dataset (`*_genomic_features.csv`)

Derived from `.gbff` files. Each row corresponds to a gene or feature entry.

| Field         | Description                                                 |
|---------------|-------------------------------------------------------------|
| Contig        | ID of the chromosome/contig                                |
| Feature_Type  | Type (gene, ncRNA, source, etc.)                           |
| Start/End     | Genomic coordinates                                        |
| Strand        | +1 or -1 orientation                                       |
| Locus_Tag     | Unique identifier                                          |
| Gene          | Gene name                                                  |
| Product       | RNA/protein description                                    |
| Protein_ID    | Protein identifier                                         |
| Translation   | Amino acid sequence (if applicable)                        |
| Note          | Comments or method used                                    |

File size: **35–75 MB**

---

##### b) Protein FASTA Dataset (`*_proteins.fasta`)

Derived from the genomic features CSV. Used in downstream tools like InterProScan.

- FASTA format (`>XP_XXXXXX...` headers with protein sequence)
- File size: **18–46 MB**

---

##### c) InterProScan Annotation Dataset (`*_interpro.tsv`)

Generated from InterProScan runs on protein sequences. Includes functional domains and pathway annotations.

Example:
```
XP_052517070.1 ... SSF144270 Eferin C-domain ... IPR037245 ... Reactome:R-HSA-432040
```

- File size: **3.75–10 GB**

---

##### d) Genomic CDS Features Dataset (`*_CDS.csv`)

Subset of genomic features containing only `CDS`, `gene`, and `mRNA` entries.

- Extracted from the original genomic features CSV.
- File size: **22–57 MB**

---

##### e) SuperMatrix Dataset (`core_orthologs_supermatrix.fasta`)

Concatenated amino acid alignments of orthologous proteins.

| Field           | Description                                           |
|------------------|-------------------------------------------------------|
| Species ID       | FASTA headers like `>Takin`, `>Yak`, `>Buffalo`       |
| Amino Acid Seq   | MAFFT-aligned, concatenated ortholog sequences        |
| Sequence Length  | Total length of concatenated orthologs                |
| Alignment Gaps   | Represented by `-` for alignment                      |
---

### 🛠️ Methodology

#### 1. **Genome-Wide Similarity with Mash**
- Fast estimation of genetic distance between species using k-mer sketches
- Output: `highres_distances.tsv`
- Yak & Buffalo are genetically closer than Takin

#### 2. **Protein Domain Analysis with InterProScan**
- Protein sequences translated from `.gbff`
- InterProScan run to detect domains, motifs, and GO terms
- Output: `.tsv` with domain annotations for each species

#### 3. **Functional Enrichment (GO)**
- Extracted biological processes and molecular functions per species
- Identified unique and shared gene functions

#### 4. **Ortholog Detection with ProteinOrtho**
- All-vs-all comparison of proteomes
- Output: Ortholog clusters shared across species
- Helps detect species-specific vs. core gene families

#### 5. **Gene Grouping and Product Distribution**
- Grouped genes by function (e.g., immune, signaling, cytoskeleton)
- Compared counts across species to detect expansion/loss trends

#### 6. **Phylogenetic Tree Construction**
- Built from MAFFT-aligned orthologous proteins
- Confirmed evolutionary distance (Takin most divergent)

---

### 📈 Key Findings

| Species       | Genetic Focus                                    |
|---------------|--------------------------------------------------|
| Wild Yak      | Cytoskeletal proteins, RNA-binding, cold response |
| Takin         | Immune expansion, structural and ECM genes        |
| Water Buffalo | Broad sensory, immune, growth & stress genes     |

#### Genetic Distances (Mash)
- Yak vs Buffalo: 97.13%
- Takin vs Buffalo: 94.80%
- Yak vs Takin: 94.56%

#### Domain Highlights
- **Yak**: PDZ, RRM, Spectrin (cold/hypoxia adaptations)
- **Takin**: Immunoglobulin, Fibronectin, ECM, GPCR
- **Buffalo**: Richest domain diversity (reproduction, immunity)

---

### 🗂️ File Structure

```
Gene_Feature_Extraction/
├── 1_genomic_feature_extraction/       # Extract CSV from GBFF
├── 2_overview_of_features/             # Plot gene feature stats
├── 3_gene_grouping/                    # Compare gene families
├── 4_protein_translation/              # Extract & convert to FASTA
├── 5_protein_extraction_and_analysis/  # Run & analyze InterProScan
├── 6_GOandKegg_Pathways/               # Enrich and cluster GO terms
├── 7_Gene_extraction/                  # Extract CDS-only features
├── 8_gene_visualization/               # Plot gene product overlaps
├── 9-ProteinOrtho-Orthologs_analysis/  # Core ortholog clustering
├── 10-Phylogenetic_Tree_ortholog/      # Build & visualize tree
```

---

### ▶️ How to Run

#### InterProScan:
```bash
interproscan.sh -i species_proteins.fasta -o output.tsv -f TSV
```

#### Gene Grouping:
```bash
python gene_grouping.py
```

#### Ortholog Detection:
```bash
proteinortho5.pl *.faa > myproject.proteinortho.tsv
```

#### Plot Heatmaps:
```bash
python extract_output_heatmap.py
```

---

### 🔧 Tools Used

- **Python**: Data processing and visualization
- **WSL / Linux**: Running heavy tools like InterProScan, Mash
- **R**: Functional clustering, GO enrichment
- **ProteinOrtho**, **MAFFT**, **IQTree**: For phylogeny and orthologs


