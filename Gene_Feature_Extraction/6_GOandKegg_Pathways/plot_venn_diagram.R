library(VennDiagram)
library(grid)
library(dplyr)
library(readr)

plot_go_venn_diagram <- function(csv_path, output_path = "compareCluster_GO/GO_shared_venn.png") {
  # Read the compareCluster 
  cc_df <- read_csv(csv_path, show_col_types = FALSE)
  
  if (!all(c("Cluster", "ID") %in% colnames(cc_df))) {
    stop("CSV must contain 'Cluster' and 'ID' columns.")
  }
  
  # Extract GO terms by species
  go_terms_list <- split(cc_df$ID, cc_df$Cluster)
  
  # Filter to species with data
  go_terms_list <- go_terms_list[sapply(go_terms_list, length) > 0]
  
  if (length(go_terms_list) < 2) {
    stop("Need at least two species with GO terms to draw a Venn diagram.")
  }
  
  # Generate Venn diagram
  venn.plot <- venn.diagram(
    x = go_terms_list,
    category.names = names(go_terms_list),
    filename = NULL,  # Draw on screen
    col = "black",
    fill = c("cornflowerblue", "green", "yellow")[1:length(go_terms_list)],
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.5,
    cat.fontface = "bold"
  )
  
  # Save to file
  png(filename = output_path, width = 800, height = 800)
  grid.draw(venn.plot)
  dev.off()
  
  message("✅ Venn diagram saved to: ", output_path)
}


csv_file <- "D:/Documents/Python Stuff - Programming/AMOD Big Data research project/Wild-Yak--Takin--and-High-Altitude-Bovids---Genomic-and-Geographic-Adaptations - NONGITHUB/Gene_Feature_Extraction/6_GOandKegg_Pathways/compareCluster_GO/GO_compareCluster.csv"

# Call the function
plot_go_venn_diagram(csv_file)
