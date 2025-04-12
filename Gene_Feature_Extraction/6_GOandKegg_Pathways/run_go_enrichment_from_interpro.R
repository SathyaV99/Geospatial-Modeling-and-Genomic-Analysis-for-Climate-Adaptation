# Load required libraries
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)

# function for compareCluster
run_comparecluster_enrichment <- function(file_paths, species_names, output_dir = "compareCluster_GO") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  all_mappings <- data.frame()
  
  for (i in seq_along(file_paths)) {
    df <- read.delim(file_paths[i], header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    colnames(df)[1:15] <- c("Protein", "MD5", "Length", "Analysis", "SignatureAcc", "SignatureDesc",
                            "Start", "End", "Score", "Status", "Date", "InterPro_ID",
                            "InterPro_desc", "GO", "Pathways")
    df <- df %>% 
      select(Protein, GO) %>%
      filter(!is.na(GO)) %>%
      separate_rows(GO, sep = "\\|") %>%
      mutate(Species = species_names[i])
    
    all_mappings <- rbind(all_mappings, df)
  }
  
  if (nrow(all_mappings) == 0) {
    message("X No GO terms found in any species.")
    return(NULL)
  }
  
  # Create a list of protein IDs per species
  gene_lists <- split(all_mappings$Protein, all_mappings$Species)
  
  cc <- compareCluster(
    geneCluster = gene_lists,
    fun = "enricher",
    TERM2GENE = all_mappings %>% select(GO, Protein),
    pvalueCutoff = 0.1  # Adjust threshold as needed
  )
  
  # Save enrichment results as CSV
  result_file <- file.path(output_dir, "GO_compareCluster.csv")
  write.csv(as.data.frame(cc), result_file, row.names = FALSE)
  message(":D Saved enrichment table to: ", result_file)
  
  # Generate and save dotplot of enrichment
  plot <- dotplot(cc, showCategory = 20) + ggtitle("GO Enrichment – All Species")
  plot_file <- file.path(output_dir, "GO_compareCluster_dotplot.png")
  ggsave(plot_file, plot, width = 14, height = 8, dpi = 300)
  message(":D  Saved dotplot to: ", plot_file)
  
  print(plot)
}

# File paths ==
file_paths <- c(
  "D:/Documents/Python Stuff - Programming/AMOD Big Data research project/Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations/Gene_Feature_Extraction/4_protein_translation/wildyak_full_interpro.tsv",
  "D:/Documents/Python Stuff - Programming/AMOD Big Data research project/Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations/Gene_Feature_Extraction/4_protein_translation/takin_interpro.tsv",
  "D:/Documents/Python Stuff - Programming/AMOD Big Data research project/Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations/Gene_Feature_Extraction/4_protein_translation/buffalo_interpro.tsv"
)

# LABELS
species_names <- c("WildYak", "Takin", "WaterBuffalo")

# Run the enrichment function
run_comparecluster_enrichment(file_paths, species_names)
