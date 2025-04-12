# README – Climate Data for Species Distribution Modeling (SDM)

*Project Title:*  
Species Distribution Modeling of Wild Yak (Bos mutus) and Takin (Budorcas taxicolor) in a Changing Climate

*Data Type:*  
High-resolution gridded climate data

*Data Size:*  
23.2 GB (24,93,81,51,783 bytes)

---

## Overview

This dataset contains climate and environmental variables used for species distribution modeling (SDM) of Wild Yak and Takin across Asia, covering the years *2009 to 2024, along with future projections for the year **2050* under *SSP245* and *SSP585* climate scenarios.

The data was sourced from *TerraClimate* for historical records and *WorldClim v2.1* for future projections. All files are in *NetCDF (.nc)* format and processed to match a common spatial resolution and extent.

---

## Folder Structure


climate_data/
│
├── ppt_1990_2020/            → Monthly precipitation files (2009–2024)
├── tmin_1990_2020/           → Monthly minimum temperature files (2009–2024)
├── tmax_1990_2020/           → Monthly maximum temperature files (2009–2024)
│
├── WorldClim2050_Asia/       → Future climate layers for 2050 (SSP245 and SSP585)
│
├── elevation_resampled_to_climate.tif    → Elevation layer aligned with climate grids
├── landmask_asia.tif                     → Landmask raster for Asia


---

## Variables Included

- *ppt*: Total monthly precipitation (mm)
- *tmin*: Minimum monthly temperature (°C)
- *tmax*: Maximum monthly temperature (°C)
- *elevation*: Elevation (meters above sea level)
- *landmask*: Binary raster to mask out ocean areas

---

## Temporal Range

- Historical Data: *2009–2024* (processed annually from monthly .nc files)
- Future Projections: *2050, based on IPCC's **SSP245* and *SSP585* scenarios

---

## Usage

These datasets were used to extract environmental features at species presence and pseudo-absence locations as part of a Random Forest classification model for SDM. All rasters were processed in Python using libraries like xarray, rasterio, and GDAL.

---

## Notes

- Each year includes 12 bands (one per month) for each variable.
- All files were preprocessed to ensure alignment in resolution, projection (EPSG:4326), and extent.
- File names follow the format: TerraClimate_<variable>_<year>.nc