# Load required libraries
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define the function for compareCluster enrichment
run_comparecluster_enrichment <- function(file_paths, species_names, 
                                          output_dir = "compareCluster_GO", 
                                          pval_cutoff = 0.2,      # Increased cutoff for sparse data
                                          show_top = 20) {
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  all_mappings <- data.frame(stringsAsFactors = FALSE)
  
  # Loop over each species/file
  for (i in seq_along(file_paths)) {
    cat("\n🔬 Processing species:", species_names[i], "\n📄 File:", file_paths[i], "\n")
    
    # Read the InterProScan TSV file
    df <- read.delim(file_paths[i], header = FALSE, sep = "\t", quote = "", 
                     stringsAsFactors = FALSE)
    
    # Ensure there are at least 15 columns
    if(ncol(df) < 15){
      warning("File ", file_paths[i], " does not have at least 15 columns. Skipping...")
      next
    }
    
    # Assign column names to the first 15 columns
    colnames(df)[1:15] <- c("Protein", "MD5", "Length", "Analysis", "SignatureAcc", "SignatureDesc",
                            "Start", "End", "Score", "Status", "Date", "InterPro_ID",
                            "InterPro_desc", "GO", "Pathways")
    
    # Extract Protein and GO columns, split the GO column by pipe ("|"), and tag species
    df <- df %>% 
      select(Protein, GO) %>%
      filter(!is.na(GO)) %>%
      separate_rows(GO, sep = "\\|") %>%
      mutate(Species = species_names[i])
    
    # Print number of unique proteins with GO annotations for the species
    num_unique <- length(unique(df$Protein))
    cat("Species:", species_names[i], "- Unique proteins with GO annotations:", num_unique, "\n")
    
    all_mappings <- rbind(all_mappings, df)
  }
  
  if (nrow(all_mappings) == 0) {
    message("❌ No GO terms found in any species.")
    return(NULL)
  }
  
  # Create a list of protein IDs per species
  gene_lists <- split(all_mappings$Protein, all_mappings$Species)
  
  # Run compareCluster enrichment analysis using 'enricher' on GO terms
  cc <- compareCluster(
    geneCluster = gene_lists,
    fun = "enricher",
    TERM2GENE = all_mappings %>% select(GO, Protein),
    pvalueCutoff = pval_cutoff  # using the provided cutoff; adjust as needed
  )
  
  # Save enrichment results to CSV
  result_file <- file.path(output_dir, "GO_compareCluster.csv")
  write.csv(as.data.frame(cc), result_file, row.names = FALSE)
  message("✅ Saved enrichment table to: ", result_file)
  
  # Generate and save a dotplot of the enrichment results
  plot <- dotplot(cc, showCategory = show_top) + 
    ggtitle("GO Enrichment – All Species")
  plot_file <- file.path(output_dir, "GO_compareCluster_dotplot.png")
  ggsave(plot_file, plot, width = 14, height = 8, dpi = 300)
  message("✅ Saved dotplot to: ", plot_file)
  
  print(plot)
  return(cc)
}

# File paths (using forward slashes)
file_paths <- c(
  "D:/Documents/Python Stuff - Programming/AMOD Big Data research project/Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations - NONGITHUB/Gene_Feature_Extraction/4_protein_translation/wildyak_full_interpro.tsv",
  "D:/Documents/Python Stuff - Programming/AMOD Big Data research project/Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations - NONGITHUB/Gene_Feature_Extraction/4_protein_translation/takin_interpro.tsv",
  "D:/Documents/Python Stuff - Programming/AMOD Big Data research project/Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations - NONGITHUB/Gene_Feature_Extraction/4_protein_translation/buffalo_interpro.tsv"
)

# Species labels
species_names <- c("WildYak", "Takin", "WaterBuffalo")

# Run the compareCluster enrichment function
cc_results <- run_comparecluster_enrichment(file_paths, species_names)
