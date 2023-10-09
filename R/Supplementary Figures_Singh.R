#*******************************************************************************
#*
#*
#*            Singh et al. (2009) extracted dataset (PMID: 19821440)                                                                           
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


## Load dataset for clustering ----
#dataset <- as.data.frame(read_excel("./31_Dataset - Extracted networks/2009/2009_19821440_Singh.xlsx", na = "NA"))[, -c(4:5, 16:20)]
load("./data/singh_dataset.RData")

# Remove studies with zero sample size
dataset_singh <- subset(dataset, sample.size > 0)

# Turn treatment columns from 'double' to 'character'
dataset_singh[, 2:3] <- lapply(dataset_singh[, 2:3], as.character)

# Turn character into factor for the corresponding characteristics
dataset_singh[, 7:13] <- lapply(dataset_singh[, 7:13], as.factor)
levels(dataset_singh$prior.drugs.failed)[levels(dataset_singh$prior.drugs.failed) == "dmard"] <- "DMARDs"    

# Check if typeof has been corrected
lapply(dataset_singh, typeof)


## Distribution of characteristics per comparison ----
# Get the plots per characteristic
distr_char <- distr_characteristics(input = dataset_singh,
                                    drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                                    rename_char = 
                                      list(5:13, 
                                           c("Duration", "Disease duration", "MTX use", "RA duration", "Anti-TNF",
                                             "Prior drugs fail", "Prior TNF fail", "Comb. biologic", "Biologic naive")),
                                    label_size = 5,
                                    axis_text_size = 11,
                                    axis_x_text_angle = 0)

# Save the figures per *Quantitative* characteristic
tiff("./Figures/Figure S1_Quantit Singh.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(distr_char$Duration, 
          distr_char$`Disease duration`,
          common.legend = TRUE, 
          legend = "bottom")
dev.off()

# Save the figures per *Qualitative* characteristic
tiff("./Figures/Figure S2_Qualit Singh.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(distr_char$`MTX use`, 
          distr_char$`Anti-TNF`, 
          distr_char$`Prior drugs fail`,
          distr_char$`Comb. biologic`, 
          distr_char$`Biologic naive`,
          distr_char$`RA duration`, 
          distr_char$`Prior TNF fail`, 
          common.legend = FALSE, 
          legend = "bottom")
dev.off()

#tiff("./Figures/Figure S2a_Qualit Singh.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
#ggarrange(distr_char$`MTX use`, 
#          distr_char$`RA duration`, 
#          distr_char$`Anti-TNF`, 
#          distr_char$`Prior drugs fail`,
#          common.legend = FALSE, 
#          legend = "bottom")
#dev.off()

#tiff("./Figures/Figure S2b_Qualit Singh.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
#ggarrange(distr_char$`Prior TNF fail`, 
#          distr_char$`Comb. biologic`, 
#          distr_char$`Biologic naive`, 
#          distr_char$`Pseudostudies`,
#          common.legend = TRUE, 
#          legend = "bottom")
#dev.off()


## Distribution of missing data in trials and characteristics ----
# Get the plots 
miss_char <- miss_characteristics(input = dataset_singh,
                                  drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                                  rename_char = 
                                    list(5:13, 
                                         c("Duration", "Disease duration", "MTX use", "RA duration", "Anti-TNF",
                                           "Prior drugs fail", "Prior TNF fail", "Comb. biologic", "Biologic naive")),
                                  label_size = 5,
                                  axis_title_size = 13,
                                  axis_text_size = 11,
                                  axis_x_text_angle = 0,
                                  legend_text_size = 12,
                                  legend_title_size = 13,
                                  strip_text_size = 11)

# Save the barplots of characteristic 
#tiff("./Figures/Figure S3_Missing barplots Singh.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
miss_char$Barplot_characteristics
#dev.off()
                                  
# Save the panel of barplots per characteristic and comparison
#tiff("./Figures/Figure S4_Missing panel Singh.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
miss_char$Barplot_combined
#dev.off()

# Save the tile plot 
#tiff("./Figures/Figure S5_Missing tileplots Singh.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
miss_char$Tileplot
#dev.off()


## Conduct *heuristic* clustering of comparisons ----
singh <- comp_clustering(input = dataset_singh, 
                         drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),  
                         rename_char = 
                           list(5:13, 
                                c("Duration", "Disease duration", "MTX use", "RA duration", "Anti-TNF",
                                  "Prior drugs fail", "Prior TNF fail", "Comb. biologic", "Biologic naive")),
                         height = FALSE,
                         optimal_clusters = 3,
                         get_plots = TRUE,
                         label_size = 4,
                         title_size = 14,
                         axis_title_size = 14,
                         axis_text_size = 14,
                         axis_x_text_angle = 0,
                         legend_text_size = 13)

# Save the panel of internal measures
tiff("./Figures/Figure S3_Internal measures Singh.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
singh$Internal_measures_panel
dev.off()

# Save the silhouette plots for 2 (connectivity index) and 3 clusters (average silhouette width & Dunn index)
tiff("./Figures/Figure S4_Silhouette plots Singh.tiff", 
     height = 30, 
     width = 60, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(comp_clustering(input = dataset_singh, 
                          drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                          optimal_clusters = 2,
                          get_plots = TRUE,
                          label_size = 4,
                          title_size = 14,
                          axis_title_size = 14,
                          axis_text_size = 14,
                          legend_text_size = 13)$Silhouette_comparisons,
          comp_clustering(input = dataset_singh, 
                          drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                          optimal_clusters = 3,
                          get_plots = TRUE,
                          label_size = 4,
                          title_size = 14,
                          axis_title_size = 14,
                          axis_text_size = 14,
                          legend_text_size = 13)$Silhouette_comparisons,
          labels = c("a)", "b)"),
          common.legend = FALSE,
          legend = "bottom")
dev.off()

# Save cophenetic coefficient in an Excel 
write_xlsx(singh$Table_cophenetic_coefficient,"./Table S5_Singh.xlsx")


## Conduct *informative* clustering of comparisons ----
singh_inf <- comp_clustering(input = dataset_singh, 
                             drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),  
                             rename_char = 
                               list(5:13, 
                                    c("Duration", "Disease duration", "MTX use", "RA duration", "Anti-TNF",
                                      "Prior drugs fail", "Prior TNF fail", "Comb. biologic", "Biologic naive")),
                             height = TRUE,
                             get_plots = TRUE)


## Distribution of characteristics per cluster ----
# Get the plots per characteristic
distr_char_cluster <- distr_characteristics(input = dataset_singh,
                                            drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                                            rename_char = 
                                              list(5:13, 
                                                   c("Duration", "Disease duration", "MTX use", "RA duration", "Anti-TNF",
                                                     "Prior drugs fail", "Prior TNF fail", "Comb. biologic", "Biologic naive")),
                                            cluster = singh,
                                            label_size = 5,
                                            axis_text_size = 12,
                                            axis_x_text_angle = 0)

# Save the figures per *Quantitative* characteristic
tiff("./Figures/Figure S5_Quantit Singh.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(distr_char_cluster$Duration, 
          distr_char_cluster$`Disease duration`,
          common.legend = TRUE, 
          legend = "bottom")
dev.off()

# Save the figures per *Qualitative* characteristic
tiff("./Figures/Figure 6_Qualit Singh.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
ggarrange(distr_char_cluster$`MTX use`, 
          distr_char_cluster$`Anti-TNF`, 
          distr_char_cluster$`Prior drugs fail`,
          distr_char_cluster$`Comb. biologic`, 
          distr_char_cluster$`Biologic naive`,
          distr_char_cluster$`RA duration`, 
          distr_char_cluster$`Prior TNF fail`,
          common.legend = FALSE, 
          legend = "bottom")
dev.off()

#tiff("./Figures/Figure S9 2b_Qualit Singh.tiff", 
#     height = 30, 
#     width = 50, 
#     units = "cm", 
#     compression = "lzw", 
#     res = 300)
#ggarrange(distr_char_cluster$`Prior TNF fail`, 
#          distr_char_cluster$`Comb. biologic`, 
#          distr_char_cluster$`Biologic naive`,
#          common.legend = TRUE, 
#          legend = "bottom")
#dev.off()