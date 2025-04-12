# extract_adaptive_go_terms.R

# Load required libraries
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(VennDiagram)
library(grid)
library(readr)

# function to run compareCluster enrichment and then extract adaptive GO terms
run_comparecluster_enrichment <- function(file_paths, species_names, 
                                          output_dir = "compareCluster_GO", 
                                          enrichment_pval = 0.1, 
                                          adaptive_pval_cutoff = 0.1,
                                          show_top = 20) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  all_mappings <- data.frame(stringsAsFactors = FALSE)
  
  # Process each species file
  for (i in seq_along(file_paths)) {
    this_file <- file_paths[i]
    species <- species_names[i]
    
    if (!file.exists(this_file)) {
      message("X File not found for species ", species, ": ", this_file)
      next
    }
    
    cat("\n Processing species:", species, "\n File:", this_file, "\n")
    
    df <- tryCatch({
      read.delim(this_file, header = FALSE, sep = "\t", quote = "", 
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      message("X Error reading file for ", species, ": ", e)
      return(NULL)
    })
    
    if (is.null(df) || ncol(df) < 15) {
      warning("File for ", species, " does not have at least 15 columns. Skipping...")
      next
    }
    
    colnames(df)[1:15] <- c("Protein", "MD5", "Length", "Analysis", "SignatureAcc", "SignatureDesc",
                            "Start", "End", "Score", "Status", "Date", "InterPro_ID",
                            "InterPro_desc", "GO", "Pathways")
    
    df <- df %>% 
      select(Protein, GO) %>%
      filter(!is.na(GO)) %>%
      separate_rows(GO, sep = "\\|") %>%
      mutate(Species = species)
    
    cat("Species:", species, "- Unique proteins with GO annotations:", length(unique(df$Protein)), "\n")
    
    all_mappings <- rbind(all_mappings, df)
  }
  
  if (nrow(all_mappings) == 0) {
    message("X No GO terms found in any species. Exiting.")
    return(NULL)
  }
  
  # Create a list of protein IDs for each species
  gene_lists <- split(all_mappings$Protein, all_mappings$Species)
  
  cc <- compareCluster(
    geneCluster = gene_lists,
    fun = "enricher",
    TERM2GENE = all_mappings %>% select(GO, Protein),
    pvalueCutoff = enrichment_pval
  )
  
  # Save the compareCluster object for future use
  save(cc, file = file.path(output_dir, "cc_results.RData"))
  
  # Save the compareCluster results to CSV
  result_file <- file.path(output_dir, "GO_compareCluster.csv")
  write.csv(as.data.frame(cc), result_file, row.names = FALSE)
  message(":D Saved enrichment table to: ", result_file)
  
  # GEN AND save dotplot of the overall enrichment results
  main_plot <- dotplot(cc, showCategory = show_top) + ggtitle("GO Enrichment – All Species")
  main_plot_file <- file.path(output_dir, "GO_compareCluster_dotplot.png")
  ggsave(main_plot_file, main_plot, width = 14, height = 8, dpi = 300)
  message(":D Saved dotplot to: ", main_plot_file)
  
  print(main_plot)
  
  # ----- Adaptive GO Terms Extraction -----
  cc_res <- cc@compareClusterResult
  
  adaptive_terms <- cc_res %>%
    filter(grepl("hypoxia|HIF|oxygen|cold|heat shock", Description, ignore.case = TRUE))
  
  if (nrow(adaptive_terms) == 0) {
    message("X No adaptive GO terms were found based on the given keywords.")
  } else {
    # Save adaptive terms to CSV
    adaptive_csv <- file.path(output_dir, "adaptive_go_terms.csv")
    write.csv(adaptive_terms, adaptive_csv, row.names = FALSE)
    message(":D Saved adaptive GO terms to: ", adaptive_csv)
    
    # For the adaptive dotplot
    if (!"geneID" %in% colnames(adaptive_terms)) {
      message("X Column 'geneID' not found in adaptive terms. Skipping adaptive dotplot generation.")
    } else {
      adaptive_gene_ids <- unlist(strsplit(adaptive_terms$geneID, split = "/"))
      adaptive_enrichment <- enricher(
        gene = adaptive_gene_ids,
        TERM2GENE = adaptive_terms[, c("ID", "geneID")],
        pvalueCutoff = adaptive_pval_cutoff
      )
      if (!is.null(adaptive_enrichment) && nrow(as.data.frame(adaptive_enrichment)) > 0) {
        adaptive_dotplot <- dotplot(adaptive_enrichment, showCategory = show_top) +
          ggtitle("Adaptive GO Terms (Hypoxia/Cold/Stress)")
        adaptive_dotplot_file <- file.path(output_dir, "adaptive_go_dotplot.png")
        ggsave(adaptive_dotplot_file, adaptive_dotplot, width = 10, height = 7, dpi = 300)
        message(":D Saved adaptive dotplot to: ", adaptive_dotplot_file)
      } else {
        message("X Adaptive enrichment did not yield significant terms.")
      }
    }
    
    # Venn diagram of GO term IDs by species
    go_terms_list <- split(cc_res$ID, cc_res$Cluster)
    # Remove any empty entries
    go_terms_list <- go_terms_list[sapply(go_terms_list, length) > 0]
    
    if (length(go_terms_list) >= 2) {
      venn_plot <- venn.diagram(
        x = go_terms_list,
        category.names = names(go_terms_list),
        filename = NULL,
        output = TRUE,
        col = "black",
        fill = c("cornflowerblue", "forestgreen", "gold")[1:length(go_terms_list)],
        alpha = 0.5,
        cex = 1.5,
        cat.cex = 1.5
      )
      
      venn_file <- file.path(output_dir, "GO_shared_venn.png")
      png(filename = venn_file, width = 800, height = 800)
      grid.draw(venn_plot)
      dev.off()
      message(":D Saved Venn diagram to: ", venn_file)
    } else {
      message("X Not enough species with GO terms to generate a Venn diagram.")
    }
  }
  
  return(cc)
}

# FILE PATHS
file_paths <- c(
  "D:/Documents/Python Stuff - Programming/AMOD Big Data research project/Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations - NONGITHUB/Gene_Feature_Extraction/4_protein_translation/wildyak_full_interpro.tsv",
  "D:/Documents/Python Stuff - Programming/AMOD Big Data research project/Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations - NONGITHUB/Gene_Feature_Extraction/4_protein_translation/takin_interpro.tsv",
  "D:/Documents/Python Stuff - Programming/AMOD Big Data research project/Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations - NONGITHUB/Gene_Feature_Extraction/4_protein_translation/buffalo_interpro.tsv"
)

# Species labels
species_names <- c("WildYak", "Takin", "WaterBuffalo")

#compareCluster enrichment + GO terms
cc_results <- run_comparecluster_enrichment(file_paths, species_names)
