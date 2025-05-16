# Genomic and Geographic Adaptations to Climate Change: A Comparative Study of Wild Yak, Takin, and High-Altitude Bovids 

This project investigates how high-altitude bovids—Wild Yak, Takin, and Water Buffalo—adapt to climate change. It combines species distribution modeling and comparative genomics to find out:

* Where their habitats are now and where they'll likely shift in the future

* What genes and protein functions help them survive in extreme environments

It uses Python and R with machine learning (Random Forest), spatial analysis (centroid tracking, jittering), and genomic tools like Mash, InterProScan, ProteinOrtho, and GO enrichment.

**Key finding:**
* Wild yak habitats are shrinking and moving uphill. Takin habitats are expanding. These results help guide future conservation efforts.

![image](https://github.com/user-attachments/assets/d7a2b6ca-a024-4243-a2f0-821560fc8903)

![image](https://github.com/user-attachments/assets/668e54af-9175-413b-8cba-680e2d794254)

![image](https://github.com/user-attachments/assets/dfb5b206-fa64-4ac6-9397-ed88ebc2df1b)


* **Repository:** Integrated SDM & Comparative Genomics pipeline for high-altitude bovids

* **Tools:** Python (scikit-learn, xarray, cartopy, geopy, GeoPandas, pandas, rasterio, GDAL), R, Shell Scripts, 

* **Techniques:** Random Forest, Spatial Jittering, Centroid Tracking

* **Data:** Occurrence records, Environmental (TerraClimate, WorldClim), Genomic assemblies

* **Genomic Tools:** Mash, InterProScan, ProteinOrtho, GO Enrichment

* **Focus:** Climate Change, High-Altitude Adaptation, Conservation

**This is the actual folder/file structure for the project. Due to the file size. Many of these files have been removed!**
```text
Genomic_and_Geographic_analysis_of_high_altitude_bovids- NONGITHUB/
├── GENETIC DATA_FROM_NCBI/
│   ├── budocras_taxicolor/
│   │   ├── ncbi_dataset/
│   │   │   └── data/
│   │   │       ├── GCA_023091745.2/
│   │   │       │   ├── GCA_023091745.2_Takin1.1_genomic.fna
│   │   │       │   └── genomic.gbffw
│   │   │       ├── GCF_023091745.1/
│   │   │       │   ├── GCF_023091745.1_Takin1.1_genomic.fna
│   │   │       │   ├── cds_from_genomic.fna
│   │   │       │   ├── genomic.gbff
│   │   │       │   ├── genomic.gff
│   │   │       │   └── genomic.gtf
│   │   │       ├── assembly_data_report.jsonl
│   │   │       ├── data_summary.tsv
│   │   │       └── dataset_catalog.json
│   │   ├── README.md
│   │   └── md5sum.txt
│   ├── water_buffalo/
│   │   ├── ncbi_dataset/
│   │   │   └── data/
│   │   │       ├── GCA_019923935.1/
│   │   │       │   ├── GCA_019923935.1_NDDB_SH_1_genomic.fna
│   │   │       │   └── genomic.gbff
│   │   │       ├── GCF_019923935.1/
│   │   │       │   ├── GCF_019923935.1_NDDB_SH_1_genomic.fna
│   │   │       │   ├── cds_from_genomic.fna
│   │   │       │   ├── genomic.gbff
│   │   │       │   ├── genomic.gff
│   │   │       │   └── genomic.gtf
│   │   │       ├── assembly_data_report.jsonl
│   │   │       ├── data_summary.tsv
│   │   │       └── dataset_catalog.json
│   │   ├── README.md
│   │   └── md5sum.txt
│   ├── wildyak/
│   │   ├── ncbi_dataset/
│   │   │   └── data/
│   │   │       ├── GCA_027580195.2/
│   │   │       │   ├── GCA_027580195.2_NWIPB_WYAK_1.1_genomic.fna
│   │   │       │   └── genomic.gbff
│   │   │       ├── GCF_027580195.1/
│   │   │       │   ├── GCF_027580195.1_NWIPB_WYAK_1.1_genomic.fna
│   │   │       │   ├── cds_from_genomic.fna
│   │   │       │   ├── genomic.gbff
│   │   │       │   ├── genomic.gff
│   │   │       │   └── genomic.gtf
│   │   │       ├── assembly_data_report.jsonl
│   │   │       ├── data_summary.tsv
│   │   │       └── dataset_catalog.json
│   │   ├── README.md
│   │   └── md5sum.txt
│   ├── budocras_taxicolor.zip
│   ├── water_buffalo.zip
│   └── wildyak.zip
├── Gene_Feature_Extraction/
│   ├── 10-Phylogenetic_Tree_ortholog/
│   │   ├── aligned_fastas/
│   │   │   └── group_3_aligned.fasta
│   │   ├── group_fastas/
│   │   │   ├── group_2.fasta
│   │   │   └── group_3.fasta
│   │   ├── circular_phylo_tree.py
│   │   ├── core_orthologs_supermatrix.fasta
│   │   ├── create_supere_matrix.py
│   │   ├── my_bovid_tree.png
│   │   ├── phylo_tree.log
│   │   ├── phylo_tree.model.gz
│   │   ├── phylo_tree_noBB.bionj
│   │   ├── phylo_tree_noBB.ckp.gz
│   │   ├── phylo_tree_noBB.iqtree
│   │   ├── phylo_tree_noBB.log
│   │   ├── phylo_tree_noBB.mldist
│   │   ├── phylo_tree_noBB.model.gz
│   │   ├── phylo_tree_noBB.treefile
│   │   └── visualize_phylo_tree.py
│   ├── 1_genomic_feature_extraction/
│   │   ├── ExtractGBFFfiles.py
│   │   ├── takin_genomic_features.csv
│   │   ├── waterbuffalo_genomic_features.csv
│   │   └── wildyak_genomic_features.csv
│   ├── 2_overview_of_features/
│   │   ├── overview_features.py
│   │   ├── takin_genomic_features_cds_plot.png
│   │   ├── takin_genomic_features_summary.txt
│   │   ├── waterbuffalo_genomic_features_cds_plot.png
│   │   ├── waterbuffalo_genomic_features_summary.txt
│   │   ├── wildyak_genomic_features_cds_plot.png
│   │   └── wildyak_genomic_features_summary.txt
│   ├── 3_gene_grouping/
│   │   ├── gene_family_output/
│   │   │   ├── gene_family_comparison.csv
│   │   │   ├── gene_family_heatmap.png
│   │   │   └── sharedfunctions_gene_family_comparison.csv
│   │   ├── gene_grouping.py
│   │   ├── gene_grouping_2.py
│   │   └── tempCodeRunnerFile.python
│   ├── 4_protein_translation/
│   │   ├── protein_fastas/
│   │   │   ├── takin_proteins.fasta
│   │   │   ├── water_buffalo_proteins.fasta
│   │   │   ├── wild_yak_proteins.fasta
│   │   │   └── wild_yak_subset.fasta
│   │   ├── buffalo_interpro.tsv
│   │   ├── protein_extract_fasta.py
│   │   ├── protein_translation_to_fasta.py
│   │   ├── takin_interpro.tsv
│   │   ├── wildyak_full_interpro.tsv
│   │   └── wildyak_subset_interpro.tsv
│   ├── 5_protein_extraction_and_analysis/
│   │   ├── interpro_output/
│   │   │   ├── interpro_domain_comparison.csv
│   │   │   └── interpro_domain_heatmap.png
│   │   ├── OLD_interpro_domain_comparison.csv
│   │   ├── OLD_interpro_domain_heatmap.png
│   │   └── extract_output_heatmap.py
│   ├── 6_GOandKegg_Pathways/
│   │   ├── compareCluster_GO/
│   │   │   ├── GO_compareCluster - Copy.csv
│   │   │   ├── GO_compareCluster.csv
│   │   │   ├── GO_compareCluster_dotplot.png
│   │   │   ├── GO_shared_venn.png
│   │   │   └── cc_results.RData
│   │   ├── interpro_output/
│   │   │   ├── functional_clustering_output/
│   │   │   │   ├── cluster_mean_heatmap.png
│   │   │   │   ├── interpro_clustered.csv
│   │   │   │   └── interpro_hierarchical_clustering.png
│   │   │   ├── interpro_domain_comparison.csv
│   │   │   └── interpro_domain_heatmap.png
│   │   ├── .RData
│   │   ├── .Rhistory
│   │   ├── combine_go_uniprot_json.py
│   │   ├── extract_adaptive_go_terms.R
│   │   ├── fetch_api_concurrently_go.py
│   │   ├── func_clustering_interpro.py
│   │   ├── go_results.json
│   │   ├── merged_enriched_GO.csv
│   │   ├── plot_venn_diagram.R
│   │   ├── run_go_enrichment_from_interpro.R
│   │   └── uniprot_results.json
│   ├── 7_Gene_extraction/
│   │   ├── extract_genes_from_csv.py
│   │   ├── takin_CDS.csv
│   │   ├── waterbuffalo_CDS.csv
│   │   └── wildyak_CDS.csv
│   ├── 8_gene_visualization/
│   │   ├── gene_product_analysis_output/
│   │   │   ├── Takin/
│   │   │   │   ├── only_takin_not_buffalo.csv
│   │   │   │   └── unique_to_takin.csv
│   │   │   ├── Wildyak/
│   │   │   │   ├── Genes_unique_to_Wildyak.xlsx
│   │   │   │   ├── only_yak_not_buffalo.csv
│   │   │   │   └── unique_to_wildyak.csv
│   │   │   ├── common genes/
│   │   │   │   ├── MAIN_common_genes_related.xlsx
│   │   │   │   ├── common_gene_product_species_matrix.csv
│   │   │   │   ├── common_gene_product_species_matrix.xlsx
│   │   │   │   ├── gene_product_shared_all_species.csv
│   │   │   │   ├── gene_product_species_matrix.csv
│   │   │   │   └── shared_all_species.csv
│   │   │   ├── waterbuffalo/
│   │   │   │   ├── only_buffalo_not_yak.csv
│   │   │   │   └── unique_to_waterbuffalo.csv
│   │   │   ├── gene_common_all_sp.txt
│   │   │   ├── gene_product_distribution_barplot.png
│   │   │   ├── gene_product_distribution_barplot_with_percentage.png
│   │   │   ├── gene_product_shared_all_species.csv
│   │   │   ├── gene_product_species_heatmap.png
│   │   │   ├── gene_product_species_matrix.csv
│   │   │   ├── gene_product_summary.csv
│   │   │   ├── gene_product_venn.png
│   │   │   ├── only_buffalo_not_yak.csv
│   │   │   ├── only_takin_not_buffalo.csv
│   │   │   ├── only_yak_not_buffalo.csv
│   │   │   ├── shared_all_species.csv
│   │   │   ├── unique_to_takin.csv
│   │   │   ├── unique_to_waterbuffalo.csv
│   │   │   └── unique_to_wildyak.csv
│   │   ├── analyze_gene_product_overlap.py
│   │   ├── gene_plot_distribution.py
│   │   └── gene_product_distribution_barplot_with_percentage.png
│   ├── 9-ProteinOrtho-Orthologs_analysis/
│   │   ├── core_ortholog_output/
│   │   │   └── core_orthologs_all_species.csv
│   │   ├── gene_family_variation_output/
│   │   │   ├── gene_family_cv_histogram.png
│   │   │   ├── gene_family_stats.tsv
│   │   │   ├── high_variation_gene_families.tsv
│   │   │   └── scatter_wild_yak_vs_takin.png
│   │   ├── gene_family_variation_output1/
│   │   │   ├── gene_family_cv_histogram.png
│   │   │   ├── gene_family_stats.tsv
│   │   │   ├── high_variation_gene_families.tsv
│   │   │   └── scatter_wild_yak_vs_takin.png
│   │   ├── orthogroup_visualization/
│   │   │   ├── barplot_core_vs_accessory.png
│   │   │   ├── barplot_orthogroup_summary.png
│   │   │   ├── orthogroup_heatmap.png
│   │   │   ├── orthogroup_presence_matrix.csv
│   │   │   ├── pca_orthogroup_presence.png
│   │   │   └── venn_orthogroups.png
│   │   ├── proteinortho/
│   │   │   ├── myproject.blast-graph
│   │   │   ├── myproject.info
│   │   │   ├── myproject.proteinortho-graph
│   │   │   ├── myproject.proteinortho-graph.summary
│   │   │   ├── myproject.proteinortho.html
│   │   │   ├── myproject.proteinortho.tsv
│   │   │   ├── takin_proteins.fasta.diamond.dmnd
│   │   │   ├── takin_proteins.fasta.len
│   │   │   ├── water_buffalo_proteins.fasta.diamond.dmnd
│   │   │   ├── water_buffalo_proteins.fasta.len
│   │   │   ├── wild_yak_proteins.fasta.diamond.dmnd
│   │   │   └── wild_yak_proteins.fasta.len
│   │   ├── analyze_genefam_var.py
│   │   ├── cafe_input_gene_counts_final_fixed_final3.tsv
│   │   ├── extract_core_orthologs_from_proteinortho.py
│   │   ├── gene_presence_barplot.png
│   │   ├── gene_presence_pca.png
│   │   ├── intermediate.csv
│   │   ├── orthogroup_analyze.py
│   │   ├── visualize_orthogroups_distr.py
│   │   └── visualize_presence_matrix_pca.py
│   ├── takin_genomic.gbff
│   ├── waterbuffalo_genomic.gbff
│   └── wildyak_genomic.gbff
├── MASH_RUN/
│   ├── distances_file/
│   │   ├── desktop.ini
│   │   ├── genomic_distance_heatmap.png
│   │   ├── heatmap_genomic_distance.py
│   │   ├── highres_distances.tsv
│   │   └── tempCodeRunnerFile.py
│   ├── GCF_019923935.1_NDDB_SH_1_genomic.fna
│   ├── GCF_023091745.1_Takin1.1_genomic.fna
│   └── GCF_027580195.1_NWIPB_WYAK_1.1_genomic.fna
├── SDM/
│   ├── Code/
│   │   └── SDM_Final.ipynb
│   ├── Data/
│   │   ├── elevation_resampled_to_climate.tif
│   │   ├── evaluation_summary.xls
│   │   ├── takin_Final_cleaned.xls
│   │   └── wild_yak_Final_cleaned.xls
│   ├── Output/
│   │   ├── sdm_takin/
│   │   │   ├── .ipynb_checkpoints/
│   │   │   │   ├── takin_suitability_2009_2009-checkpoint.png
│   │   │   │   ├── takin_suitability_2010_2010-checkpoint.png
│   │   │   │   ├── takin_suitability_2018_2018-checkpoint.png
│   │   │   │   ├── takin_suitability_2020_2020-checkpoint.png
│   │   │   │   └── takin_suitability_2024_2024-checkpoint.png
│   │   │   ├── centroid_arrow_clean_map_takin.png
│   │   │   ├── centroid_distance_changes_takin.csv
│   │   │   ├── centroid_shift_trends_takin.png
│   │   │   ├── centroid_shifts_takin.csv
│   │   │   ├── landmask_asia_cropped.npy
│   │   │   ├── random_forest_model_takin.pkl
│   │   │   ├── suitability_map_2009_2009.npy
│   │   │   ├── suitability_map_2010_2010.npy
│   │   │   ├── suitability_map_2011_2011.npy
│   │   │   ├── suitability_map_2012_2012.npy
│   │   │   ├── suitability_map_2013_2013.npy
│   │   │   ├── suitability_map_2014_2014.npy
│   │   │   ├── suitability_map_2015_2015.npy
│   │   │   ├── suitability_map_2016_2016.npy
│   │   │   ├── suitability_map_2017_2017.npy
│   │   │   ├── suitability_map_2018_2018.npy
│   │   │   ├── suitability_map_2019_2019.npy
│   │   │   ├── suitability_map_2020_2020.npy
│   │   │   ├── suitability_map_2021_2021.npy
│   │   │   ├── suitability_map_2022_2022.npy
│   │   │   ├── suitability_map_2023_2023.npy
│   │   │   ├── suitability_map_2024_2024.npy
│   │   │   ├── suitability_map_2050_SSP245.npy
│   │   │   ├── suitability_map_2050_SSP585.npy
│   │   │   ├── takin_himalaya_area_trend.png
│   │   │   ├── takin_suitability_2009_2009.png
│   │   │   ├── takin_suitability_2010_2010.png
│   │   │   ├── takin_suitability_2011_2011.png
│   │   │   ├── takin_suitability_2012_2012.png
│   │   │   ├── takin_suitability_2013_2013.png
│   │   │   ├── takin_suitability_2014_2014.png
│   │   │   ├── takin_suitability_2015_2015.png
│   │   │   ├── takin_suitability_2016_2016.png
│   │   │   ├── takin_suitability_2017_2017.png
│   │   │   ├── takin_suitability_2018_2018.png
│   │   │   ├── takin_suitability_2019_2019.png
│   │   │   ├── takin_suitability_2020_2020.png
│   │   │   ├── takin_suitability_2021_2021.png
│   │   │   ├── takin_suitability_2022_2022.png
│   │   │   ├── takin_suitability_2023_2023.png
│   │   │   ├── takin_suitability_2024_2024.png
│   │   │   ├── takin_suitability_2050_SSP245.png
│   │   │   ├── takin_suitability_2050_SSP585.png
│   │   │   └── takin_suitability_area_trend.png
│   │   └── sdm_yak/
│   │       ├── .ipynb_checkpoints/
│   │       │   ├── area_trend-checkpoint.png
│   │       │   ├── centroid_distance_changes-checkpoint.csv
│   │       │   ├── centroid_lat_trend-checkpoint.png
│   │       │   ├── centroid_lon_trend-checkpoint.png
│   │       │   ├── centroid_shift_trends-checkpoint.png
│   │       │   ├── centroid_shifts-checkpoint.csv
│   │       │   ├── centroid_trends_plotly-checkpoint.html
│   │       │   ├── elevation_trend-checkpoint.png
│   │       │   ├── habitat_gain_loss_2050_SSP245-checkpoint.png
│   │       │   ├── habitat_gain_loss_2050_SSP585_fixed-checkpoint.png
│   │       │   ├── habitat_gain_loss_map_2050_SSP245-checkpoint.png
│   │       │   ├── habitat_gain_loss_map_SSP245_cartopy-checkpoint.png
│   │       │   ├── habitat_gain_loss_map_SSP585_2050_final-checkpoint.png
│   │       │   ├── yak_suitability_2009_2009-checkpoint.png
│   │       │   ├── yak_suitability_2011_2011-checkpoint.png
│   │       │   ├── yak_suitability_2012_2012-checkpoint.png
│   │       │   ├── yak_suitability_2019_2019-checkpoint.png
│   │       │   ├── yak_suitability_2020_2020-checkpoint.png
│   │       │   ├── yak_suitability_2024_2024-checkpoint.png
│   │       │   ├── yak_suitability_2050_SSP245-checkpoint.png
│   │       │   ├── yak_suitability_2050_SSP585-checkpoint.png
│   │       │   └── yak_suitability_area_trend_final-checkpoint.png
│   │       ├── centroid_arrow_clean_map-Copy1.png
│   │       ├── centroid_arrow_clean_map.png
│   │       ├── centroid_distance_changes.csv
│   │       ├── centroid_shift_trends.png
│   │       ├── centroid_shifts.csv
│   │       ├── centroid_trends_plotly.html
│   │       ├── landmask.npy
│   │       ├── landmask_asia_cropped.npy
│   │       ├── random_forest_model.pkl
│   │       ├── suitability_map_2009_2009.npy
│   │       ├── suitability_map_2010_2010.npy
│   │       ├── suitability_map_2011_2011.npy
│   │       ├── suitability_map_2012_2012.npy
│   │       ├── suitability_map_2013_2013.npy
│   │       ├── suitability_map_2014_2014.npy
│   │       ├── suitability_map_2015_2015.npy
│   │       ├── suitability_map_2016_2016.npy
│   │       ├── suitability_map_2017_2017.npy
│   │       ├── suitability_map_2018_2018.npy
│   │       ├── suitability_map_2019_2019.npy
│   │       ├── suitability_map_2020_2020.npy
│   │       ├── suitability_map_2021_2021.npy
│   │       ├── suitability_map_2022_2022.npy
│   │       ├── suitability_map_2023_2023.npy
│   │       ├── suitability_map_2024_2024.npy
│   │       ├── suitability_map_2050_SSP245.npy
│   │       ├── suitability_map_2050_SSP585.npy
│   │       ├── yak_suitability_2009_2009-Copy1.png
│   │       ├── yak_suitability_2009_2009.png
│   │       ├── yak_suitability_2010_2010.png
│   │       ├── yak_suitability_2011_2011.png
│   │       ├── yak_suitability_2012_2012.png
│   │       ├── yak_suitability_2013_2013.png
│   │       ├── yak_suitability_2014_2014.png
│   │       ├── yak_suitability_2015_2015.png
│   │       ├── yak_suitability_2016_2016.png
│   │       ├── yak_suitability_2017_2017.png
│   │       ├── yak_suitability_2018_2018.png
│   │       ├── yak_suitability_2019_2019.png
│   │       ├── yak_suitability_2020_2020.png
│   │       ├── yak_suitability_2021_2021.png
│   │       ├── yak_suitability_2022_2022.png
│   │       ├── yak_suitability_2023_2023.png
│   │       ├── yak_suitability_2024_2024.png
│   │       ├── yak_suitability_2050_SSP245-Copy1.png
│   │       ├── yak_suitability_2050_SSP245.png
│   │       ├── yak_suitability_2050_SSP585.png
│   │       └── yak_suitability_area_trend_final.png
│   └── README.md
├── .gitattributes
├── folder_tree.py
├── ~$README.docx
└── ~WRL0830.tmp
```
