#*******************************************************************************
#*
#*
#*            Baker et al. (2009) extracted dataset (PMID: 19637942)                                                                           
#*
#*
#*******************************************************************************


## Load rnmamod developmental version
#remotes::install_github("https://github.com/LoukiaSpin/rnmamod.git", force = TRUE)


## Load libraries 
list.of.packages <- c("readxl", "ggplot2", "reshape2", "ggpubr", "heatmaply", 
                      "igraph", "gridExtra", "writexl", "scales", "stringr", "rnmamod")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)


## Load functions ----
# Plot distribution of characteristics
source("./R/distr.characteristics_function.R")

# Visualise missingness in characteristics
source("./R/miss.characteristics_function.R")

# Gower dissimilarity
source("./R/gower.distance_function.R")

# Cluster comparisons
source("./R/comp.clustering_function.R")

# Connectivity index
source("./R/connectivity.index_function.R")

# Silhouette index
source("./R/silhouette.index_function.R")

# Dunn index
source("./R/dunn.index_function.R")

# Panel of all internal measures 
source("./R/internal.measures.plot_function.R")

# Network of comparisons 
source("./R/network.comparisons_function.R")

# Dendrogram with heatmap
source("./R/dendro.heatmap_function.R")

# Turn long (netmeta) into wide (gemtc) format
source("./R/long.to.wide_function.R")


## Load dataset for clustering ----
#dataset <- as.data.frame(read_excel("./31_Dataset - Extracted networks/2009/2009_19637942_Baker.xlsx", na = "NA"))[, -c(4:5, 15, 18:21)]
load("./data/baker_dataset.RData")

# Turn treatment columns from 'double' to 'character'
dataset[, 2:3] <- lapply(dataset[, 2:3], as.character)

# Turn character into factor for the corresponding characteristics
dataset[, 6:8] <- lapply(dataset[, 6:8], as.factor)

# Check if typeof has been corrected
lapply(dataset, typeof)


## Load dataset for network plot ----
#data_network <- as.data.frame(read_excel("./31_Dataset - Extracted networks/2009/2009_19637942_Baker.xlsx", na = "NA"))[, c(1:3, 20:21)]
load("./data/baker_network.RData")

# Turn treatment columns from 'double' to 'character'
data_network[, 2:3] <- lapply(data_network[, 2:3], as.character)

# Check if typeof has been corrected
lapply(data_network, typeof)

# Turn into the proper format 
network_data <- long_to_wide(input = data_network)


## Distribution of characteristics per comparison ----
# Get the plots per characteristic
distr_char <- distr_characteristics(input = dataset,
                                    drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                                    rename_char = list(c(5:14),
                                                       c("Duration", "Random allocation", "Double blinding",
                                                         "Withdrawals description", "Quality score", "Inclusion FEV1",
                                                         "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                                    label_size = 3.0,
                                    axis_text_size = 11,
                                    axis_x_text_angle = 45)

# Save the figures per *Quantitative* characteristic 
tiff("./Figures/Figure S7_Quantit Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(distr_char$`Duration`, 
          distr_char$`Quality score`, 
          distr_char$`Inclusion FEV1`, 
          distr_char$`Inclusion FVC`,
          distr_char$`Smoking history`, 
          distr_char$`Minimum FEV1 (%)`, 
          distr_char$`Maximum FEV1 (%)`, 
          common.legend = TRUE, 
          legend = "bottom")
dev.off()

#tiff("./Figures/Figure S7a_Quantit Baker.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
#ggarrange(distr_char$`Duration`, 
#          distr_char$`Quality score`, 
#          distr_char$`Inclusion FEV1`, 
#          distr_char$`Inclusion FVC`,
#          common.legend = TRUE, 
#          legend = "bottom")
#dev.off()

# Save the figures per *Quantitative* characteristic 
#tiff("./Figures/Figure S7b_Quantit Baker.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
#ggarrange(distr_char$`Smoking history`, 
#          distr_char$`Minimum FEV1 (%)`, 
#          distr_char$`Maximum FEV1 (%)`, 
#          common.legend = TRUE, 
#          legend = "bottom")
#dev.off()

# Save the figures per *Qualitative* characteristic
tiff("./Figures/Figure S8_Qualit Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(distr_char$`Random allocation`, 
          distr_char$`Double blinding`, 
          distr_char$`Withdrawals description`, 
          distr_char$Pseudostudies,
          common.legend = FALSE, 
          legend = "bottom")
dev.off()


## Distribution of missing data in trials and characteristics ----
# Get the plots 
miss_char <- miss_characteristics(input = dataset,
                                  drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                                  rename_char = list(c(5:14),
                                                     c("Duration", "Random allocation", "Double blinding",
                                                       "Withdrawals description", "Quality score", "Inclusion FEV1",
                                                       "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                                  label_size = 4.0,
                                  axis_title_size = 11, #14,
                                  axis_text_size = 9.0, #12,
                                  axis_x_text_angle = 0,
                                  legend_text_size = 12, #13,
                                  legend_title_size = 11, #14,
                                  strip_text_size = 8.0,
                                  strip_text_angle = 0)

# Save the barplots of characteristic 
tiff("./Figures/Figure S9_Missing barplots Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
miss_char$Barplot_characteristics
dev.off()

# Save the panel of barplots per characteristic and comparison
tiff("./Figures/Figure S10_Missing panel Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
miss_char$Barplot_combined
dev.off()

# Save the tile plot 
tiff("./Figures/Figure S11_Missing tileplots Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
miss_char$Tileplot
dev.off()


## Conduct *heuristic* clustering of comparisons ----
baker <- comp_clustering(input = dataset, 
                         drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                         rename_char = list(c(5:14),
                                            c("Duration", "Random allocation", "Double blinding",
                                              "Withdrawals description", "Quality score", "Inclusion FEV1",
                                              "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                         num_neighb = 2, 
                         optimal_clusters = 5,
                         get_plots = TRUE,
                         label_size = 4.0,
                         axis_text_size = 11)

# Save the barplot on total dissimilarities
tiff("./Figures/Figure S12_Total diss barplot Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
baker$Total_dissimilarity_plot
dev.off()

# Save the panel of internal measures
tiff("./Figures/Figure S13_Internal measures Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
baker$Internal_measures_panel
dev.off()

# Get the panel of four silhouette plots on 2, 3, 4, and 5 clusters
silh2 <- comp_clustering(input = dataset, 
                         drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                         optimal_clusters = 2,
                         get_plots = TRUE,
                         label_size = 3,
                         title_size = 11,
                         axis_title_size = 11,
                         axis_text_size = 11, 
                         legend_text_size = 11)$Silhouette_comparisons
silh4 <- comp_clustering(input = dataset, 
                         drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                         optimal_clusters = 3,
                         get_plots = TRUE,
                         label_size = 3,
                         title_size = 11,
                         axis_title_size = 11,
                         axis_text_size = 11, 
                         legend_text_size = 11)$Silhouette_comparisons
silh8 <- comp_clustering(input = dataset, 
                         drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                         optimal_clusters = 4,
                         get_plots = TRUE,
                         label_size = 3,
                         title_size = 11,
                         axis_title_size = 11,
                         axis_text_size = 11,
                         legend_text_size = 11)$Silhouette_comparisons
silh11 <- comp_clustering(input = dataset, 
                          drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                          optimal_clusters = 5,
                          get_plots = TRUE,
                          label_size = 3,
                          title_size = 11,
                          axis_title_size = 11,
                          axis_text_size = 11,
                          legend_text_size = 11)$Silhouette_comparisons

# Save the panel of four silhouette plots on 2, 3, 4, and 5 clusters
tiff("./Figures/Figure S14_Silhouette plots Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(silh2, silh4, silh8, silh11)
dev.off()

# Save cophenetic coefficient in an Excel 
write_xlsx(baker$Table_cophenetic_coefficient,"./Table S6_Baker.xlsx")


## Conduct *informative* clustering of comparisons ----
baker_inf <- comp_clustering(input = dataset, 
                             drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                             rename_char = list(c(5:14),
                                                c("Duration", "Random allocation", "Double blinding",
                                                  "Withdrawals description", "Quality score", "Inclusion FEV1",
                                                  "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                             height = TRUE,
                             get_plots = TRUE,
                             label_size = 4.0,
                             axis_text_size = 11)

## Distribution of characteristics per cluster (INFORMATIVE CLUSTERING) ----
# Get the plots per characteristic
distr_char_cluster_inf <- 
  distr_characteristics(input = dataset,
                        drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                        rename_char = list(c(5:14),
                                           c("Duration", "Random allocation", "Double blinding",
                                             "Withdrawals description", "Quality score", "Inclusion FEV1",
                                             "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                        cluster = baker_inf,
                        label_size = 3.0,
                         axis_text_size = 12)

# Save the figures per *Quantitative* characteristic 
tiff("./Figures/Figure S15_Quantit Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(distr_char_cluster_inf$`Duration`, 
          distr_char_cluster_inf$`Quality score`, 
          distr_char_cluster_inf$`Inclusion FEV1`, 
          distr_char_cluster_inf$`Inclusion FVC`,
          distr_char_cluster_inf$`Smoking history`, 
          distr_char_cluster_inf$`Minimum FEV1 (%)`, 
          distr_char_cluster_inf$`Maximum FEV1 (%)`, 
          common.legend = TRUE, 
          legend = "bottom")
dev.off()

#tiff("./Figures/Figure S15_Quantit Baker.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
#ggarrange(distr_char_cluster_inf$`Duration`, 
#          distr_char_cluster_inf$`Quality score`, 
#          distr_char_cluster_inf$`Inclusion FEV1`, 
#          distr_char_cluster_inf$`Inclusion FVC`, 
#          common.legend = TRUE, 
#          legend = "bottom")
#dev.off()

# Save the figures per *Quantitative* characteristic 
#tiff("./Figures/Figure S18b_Quantit Baker.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
#ggarrange(distr_char_cluster_inf$`Smoking history`, 
#          distr_char_cluster_inf$`Minimum FEV1 (%)`, 
#          distr_char_cluster_inf$`Maximum FEV1 (%)`, 
#          common.legend = TRUE, 
#          legend = "bottom")
#dev.off()

# Save the figures per *Qualitative* characteristic
tiff("./Figures/Figure S16_Qualit Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(distr_char_cluster_inf$`Random allocation`, 
          distr_char_cluster_inf$`Double blinding`, 
          distr_char_cluster_inf$`Withdrawals description`, 
          common.legend = FALSE, 
          legend = "bottom")
dev.off()


## Distribution of characteristics per cluster (HEURISTIC CLUSTERING) ----
# Get the plots per characteristic
distr_char_cluster <- distr_characteristics(input = dataset,
                                            drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                                            rename_char = list(c(5:14),
                                                               c("Duration", "Random allocation", "Double blinding",
                                                                 "Withdrawals description", "Quality score", "Inclusion FEV1",
                                                                 "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                                            cluster = baker,
                                            label_size = 3.0,
                                            axis_text_size = 12)

# Save the figures per *Quantitative* characteristic 
tiff("./Figures/Figure S17_Quantit Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(distr_char_cluster$`Duration`, 
          distr_char_cluster$`Quality score`, 
          distr_char_cluster$`Inclusion FEV1`, 
          distr_char_cluster$`Inclusion FVC`,
          distr_char_cluster$`Smoking history`, 
          distr_char_cluster$`Minimum FEV1 (%)`, 
          distr_char_cluster$`Maximum FEV1 (%)`, 
          common.legend = TRUE, 
          legend = "bottom")
dev.off()

#tiff("./Figures/Figure S17a_Quantit Baker.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
#ggarrange(distr_char_cluster$`Duration`, 
#          distr_char_cluster$`Quality score`, 
#          distr_char_cluster$`Inclusion FEV1`, 
#          distr_char_cluster$`Inclusion FVC`,
#          common.legend = TRUE, 
#          legend = "bottom")
#dev.off()

# Save the figures per *Quantitative* characteristic 
#tiff("./Figures/Figure S17b_Quantit Baker.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
#ggarrange(distr_char_cluster$`Smoking history`, 
#          distr_char_cluster$`Minimum FEV1 (%)`, 
#          distr_char_cluster$`Maximum FEV1 (%)`, 
#          common.legend = TRUE, 
#          legend = "bottom")
#dev.off()

# Save the figures per *Qualitative* characteristic
tiff("./Figures/Figure S18_Qualit Baker.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(distr_char_cluster$`Random allocation`, 
          distr_char_cluster$`Double blinding`, 
          distr_char_cluster$`Withdrawals description`, 
          common.legend = FALSE, 
          legend = "bottom")
dev.off()


## Figure S22: Dendrogram under HEURISTIC clustering ----
# Save as png, bring to pptx and save as tiff
#dendro_heatmap(input = baker_inf)


## Figure S23: Network of comparisons under HEURISTIC clustering ----
# First, check dimensions
#layout(matrix(c(1, 2, 3, 4, 5, 5), ncol = 2, byrow = TRUE), 
#       heights = c(0.9, 0.9, 0.1),  
#       widths = c(0.9, 0.9))
#layout.show(n = 5)
#dev.off()

# Then get the plot
#tiff("./Figures/Figure S23.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
#layout(matrix(c(1, 2, 3, 4, 5, 5), ncol =  2, byrow = TRUE), heights = c(1.1, 1.1, 0.2))
#par(mar = c(0, 0, 0, 0), oma = c(1.2, 1.2, 1.2, 1.2))

#network_comparisons(data = dataset, 
#                    drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
#                    total_diss = baker$Total_dissimilarity, 
#                    comp_diss = baker$Dissimilarity_table,
#                    optimal_dist = "euclidean",
#                    optimal_link = "average",
#                    optimal_clusters = 5,
#                    node_frame_width = 1.5,
#                    node_label_cex = 1.5,
#                    node_label_dist = 0,
#                    edge_level = "low")
#mtext(expression(paste(bold("a)"))), side = 3, line = -1.6, adj = 0.20)

#network_comparisons(data = dataset, 
#                    drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
#                    total_diss = baker$Total_dissimilarity, 
#                    comp_diss = baker$Dissimilarity_table,
#                    optimal_dist = "euclidean",
#                    optimal_link = "average",
#                    optimal_clusters = 5,
#                    node_frame_width = 1.5,
#                    node_label_cex = 1.5,
#                    node_label_dist = 0,
#                    edge_level = "moderate")
#mtext(expression(paste(bold("b)"))), side = 3, line = -1.6,  adj = 0.20)

#network_comparisons(data = dataset, 
#                    drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
#                    total_diss = baker$Total_dissimilarity, 
#                    comp_diss = baker$Dissimilarity_table,
#                    optimal_dist = "euclidean",
#                    optimal_link = "average",
#                    optimal_clusters = 5,
#                    node_frame_width = 1.5,
#                    node_label_cex = 1.5,
#                    node_label_dist = 0,
#                    edge_level = "high")
#mtext(expression(paste(bold("c)"))), side = 3, line = -1.6,  adj = 0.20)

#network_comparisons(data = dataset, 
#                    drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
#                    total_diss = baker$Total_dissimilarity, 
#                    comp_diss = baker$Dissimilarity_table,
#                    optimal_dist = "euclidean",
#                    optimal_link = "average",
#                    optimal_clusters = 5,
#                    node_frame_width = 1.5,
#                    node_label_cex = 1.5,
#                    node_label_dist = 0,
#                    edge_level = "very high")
#mtext(expression(paste(bold("d)"))), side = 3, line = -1.6,  adj = 0.20)

#par(mar = c(0, 0, 0, 0)) # Make the margins 0
#plot(1, type = "n", axes = F, xlab = "", ylab = "") 

#legend("bottom", 
#       legend = c("Low [0.00, 0.25]", "Moderate (0.25, 0.50]", "High (0.50, 0.75]", "Very high (0.75, 1.00]"), 
#       col = c("#A6D854","#E6AB02", "#D95F02", "#E31A1C"), 
#       pch = 15,
#       text.width = c(0.064, 0.083, 0.065, 0.07),
#       xjust = 0, 
#       yjust = 0,
#       inset = 0.2,
#       xpd = TRUE,
#       cex = 1.7,
#       bty = "n",
#       horiz = TRUE)
#dev.off()
