#*******************************************************************************
#*
#*  R Code to replicate the Supplementary Figures & Tables 
#*  <Baker and colleagues (2009; PMID: 19637942) - COPD network>
#*     
#*  Authors: Loukia M. Spineli, Katerina Papadimitropoulou, Chrysostomos Kalyvas
#*  Date: December 2023
#*  
#*******************************************************************************


## Load rnmamod developmental version
remotes::install_github("https://github.com/LoukiaSpin/rnmamod.git", force = TRUE)


## Load libraries ----
list.of.packages <- c("rnmamod", "reshape2", "ggpubr", "stringr", "cluster", "scales", "dendextend")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)


## Load Baker dataset for clustering ----
load("./data/baker_dataset.RData")

# Turn treatment columns from 'double' to 'character'
dataset[, 2:3] <- lapply(dataset[, 2:3], as.character)

# Turn character into factor for the corresponding characteristics
dataset[, 6:8] <- lapply(dataset[, 6:8], as.factor)

# Check if typeof has been corrected
lapply(dataset, typeof)


## Distribution of characteristics per comparison ----
# ?distr_characteristics
distr_char <- distr_characteristics(input = dataset,
                                    drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                                    rename_char = list(c(5:14),
                                                       c("Duration", "Random allocation", "Double blinding",
                                                         "Withdrawals description", "Quality score", "Inclusion FEV1",
                                                         "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                                    label_size = 2.0,
                                    title_size = 11,
                                    axis_text_size = 10,
                                    axis_x_text_angle = 45)

# Figure S10
tiff("./Figure S10.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
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

# Figure S11
tiff("./Figure S11.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
ggarrange(distr_char$`Random allocation`, 
          distr_char$`Double blinding`, 
          distr_char$`Withdrawals description`, 
          distr_char$Pseudostudies,
          common.legend = FALSE, 
          legend = "bottom")
dev.off()


## Distribution of missing data in trials and characteristics ----
# ?miss_characteristics 
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

# Figure S12 
miss_char$Barplot_characteristics

# Figure S13
miss_char$Barplot_combined

# Figure S14
miss_char$Tileplot


## Create dissimilarity matrix (D) and heatmap ----
# ?comp_clustering
baker <- comp_clustering(input = dataset, 
                         weight = c(rep(1, 9), 0.5, 0.5),
                         drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                         threshold = 0.13,
                         informative = TRUE,
                         get_plots = TRUE,
                         label_size = 4.0,
                         axis_text_size = 11)


## Conduct hierarchical clustering ----
# ?comp_clustering
baker_hie <- comp_clustering(input = dataset, 
                             weight = c(rep(1, 9), 0.5, 0.5),
                             drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                             informative = FALSE,
                             optimal_clusters = 2,
                             get_plots = TRUE,
                             label_size = 4.0,
                             axis_text_size = 11)

# Table S5: Cophenetic correlation coefficients 
baker_hie$Table_cophenetic_coefficient
 
# Figure S15
tiff("./Figure S15.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
baker$Within_comparison_dissimilarity
dev.off()

# Figure S16
tiff("./Figure S16.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
baker$Between_comparison_dissimilarity
dev.off()

# Figure S17
tiff("./Figure S17.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
baker_hie$Profile_plot
dev.off()

# Figure S18a
tiff("./Figure S18a.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
comp_clustering(input = dataset, 
                weight = c(rep(1, 9), 0.5, 0.5),
                drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"), 
                informative = FALSE,
                optimal_clusters = 2,
                get_plots = TRUE,
                label_size = 2.5,
                title_size = 10,
                axis_title_size = 10,
                axis_text_size = 7,
                legend_text_size = 10)$Silhouette_width_plot
dev.off()

# Figure S18b
tiff("./Figure S18b.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
comp_clustering(input = dataset, 
                weight = c(rep(1, 9), 0.5, 0.5),
                drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"), 
                informative = FALSE,
                optimal_clusters = 3,
                get_plots = TRUE,
                label_size = 2.5,
                title_size = 10,
                axis_title_size = 10,
                axis_text_size = 7,
                legend_text_size = 10)$Silhouette_width_plot
dev.off()

# Figure S18c
tiff("./Figure S18c.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
comp_clustering(input = dataset, 
                weight = c(rep(1, 9), 0.5, 0.5),
                drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"), 
                informative = FALSE,
                optimal_clusters = 4,
                get_plots = TRUE,
                label_size = 2.5,
                title_size = 10,
                axis_title_size = 10,
                axis_text_size = 7,
                legend_text_size = 10)$Silhouette_width_plot
dev.off()

# Figure S18d
tiff("./Figure S18d.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
comp_clustering(input = dataset, 
                weight = c(rep(1, 9), 0.5, 0.5),
                drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"), 
                informative = FALSE,
                optimal_clusters = 5,
                get_plots = TRUE,
                label_size = 2.5,
                title_size = 10,
                axis_title_size = 10,
                axis_text_size = 7,
                legend_text_size = 10)$Silhouette_width_plot
dev.off()

# Figure S19
dendro_heatmap(input = baker_hie,
               label_size = 4,
               axis_text_size = 6.0)


## Distribution of characteristics per cluster ----
distr_char_cluster <- distr_characteristics(input = dataset,
                                            drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                                            rename_char = list(c(5:14),
                                                               c("Duration", "Random allocation", "Double blinding",
                                                                 "Withdrawals description", "Quality score", "Inclusion FEV1",
                                                                 "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                                            cluster = baker_hie,
                                            label_size = 4.0,
                                            axis_text_size = 12)

# Figure S20
tiff("./Figure S20.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
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

# Figure S21
tiff("./Figure S21.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
ggarrange(distr_char_cluster$`Random allocation`, 
          distr_char_cluster$`Double blinding`, 
          distr_char_cluster$`Withdrawals description`, 
          common.legend = FALSE, 
          legend = "bottom")
dev.off()
