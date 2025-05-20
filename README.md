# Genomic and Geographic Adaptations to Climate Change: A Comparative Study of Wild Yak, Takin, and High-Altitude Bovids 

## Wild Yak faces substantial habitat loss and altitudinal displacement by 2050, driven by climate change and limited adaptive genomic traits.

### Computationally intensive analyses like InterProScan and BLAST averaged 2–4 days each for a single species to complete!

This project investigates how high-altitude bovids-Wild Yak, Takin, and Water Buffalo-adapt to climate change. It combines species distribution modeling and comparative genomics to find out:

* Where their habitats are now and where they'll likely shift in the future

* What genes and protein functions help them survive in extreme environments

It uses Python and R with machine learning (Random Forest), spatial analysis (centroid tracking, jittering), and genomic tools like Mash, InterProScan, ProteinOrtho, and GO enrichment.

![image](https://github.com/user-attachments/assets/dfb5b206-fa64-4ac6-9397-ed88ebc2df1b)

**Key finding:**
* Wild yak habitats are shrinking and moving uphill. Takin habitats are expanding. These results help guide future conservation efforts.

![image](https://github.com/user-attachments/assets/668e54af-9175-413b-8cba-680e2d794254)

Additionally, wild yaks have fewer heat shock genes, and due to their thick fur, they struggle to regulate body heat. They are adapted to cold climates at elevations around 3,000 feet and cannot tolerate warmer temperatures.

<img width="1012" alt="immune_genes" src="https://github.com/user-attachments/assets/64fc52cd-2be0-4b08-83bd-89689bc10c3c" />

## 🌍 Species Distribution Modeling (SDM)

This part of the project models current and future habitat suitability for **Wild Yak** and **Takin** using geospatial and climate data. It applies machine learning to predict where these animals can survive based on environmental conditions.

---

### 🔸 Objective

Predict species range shifts from **2009 to 2050** using environmental variables and occurrence records.

---

### 📦 Data Sources

| Data Type           | Source            | Years Covered   | Format     |
|---------------------|-------------------|------------------|------------|
| Occurrence Data     | GBIF, iNaturalist, Literature | 2014–2024 | `.csv` |
| Climate Data (Current) | TerraClimate       | 2009–2024       | `.nc`      |
| Climate Projections | WorldClim (SSP245 & SSP585) | 2050 | `.nc` |
| Elevation Data      | Google Earth Engine | -                | `.tif`     |
| Landmask            | Natural Earth        | -                | `.shp` → `.tif` |

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


