#*******************************************************************************
#*
#*  R Code to replicate the Supplementary Figures & Tables 
#*  <Baker and colleagues (2009; PMID: 19637942) - COPD network>
#*     
#*  Authors: Loukia M. Spineli, Katerina Papadimitropoulou, Chrysostomos Kalyvas
#*  Date: October 2023
#*  
#*******************************************************************************


## Load rnmamod developmental version
remotes::install_github("https://github.com/LoukiaSpin/rnmamod.git", force = TRUE)


## Load libraries ----
list.of.packages <- c("rnmamod", "reshape2", "ggpubr", "stringr")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)


## Load functions ----
source("./R/long.to.wide_function.R")


## Load Baker dataset for clustering ----
load("./data/baker_dataset.RData")

# Turn treatment columns from 'double' to 'character'
dataset[, 2:3] <- lapply(dataset[, 2:3], as.character)

# Turn character into factor for the corresponding characteristics
dataset[, 6:8] <- lapply(dataset[, 6:8], as.factor)

# Check if typeof has been corrected
lapply(dataset, typeof)


## Load dataset for network plot ----
load("./data/baker_network.RData")

# Turn treatment columns from 'double' to 'character'
data_network[, 2:3] <- lapply(data_network[, 2:3], as.character)

# Check if typeof has been corrected
lapply(data_network, typeof)

# Turn into the proper format 
network_data <- long_to_wide(input = data_network)


## Distribution of characteristics per comparison ----
# ?distr_characteristics
distr_char <- distr_characteristics(input = dataset,
                                    drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                                    rename_char = list(c(5:14),
                                                       c("Duration", "Random allocation", "Double blinding",
                                                         "Withdrawals description", "Quality score", "Inclusion FEV1",
                                                         "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                                    label_size = 3.0,
                                    axis_text_size = 11,
                                    axis_x_text_angle = 45)

# Figure S7
ggarrange(distr_char$`Duration`, 
          distr_char$`Quality score`, 
          distr_char$`Inclusion FEV1`, 
          distr_char$`Inclusion FVC`,
          distr_char$`Smoking history`, 
          distr_char$`Minimum FEV1 (%)`, 
          distr_char$`Maximum FEV1 (%)`, 
          common.legend = TRUE, 
          legend = "bottom")

# Figure S8
ggarrange(distr_char$`Random allocation`, 
          distr_char$`Double blinding`, 
          distr_char$`Withdrawals description`, 
          distr_char$Pseudostudies,
          common.legend = FALSE, 
          legend = "bottom")


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

# Figure S9 
miss_char$Barplot_characteristics

# Figure S10
miss_char$Barplot_combined

# Figure S11
miss_char$Tileplot


## Conduct *HEURISTIC* clustering of comparisons ----
# ?comp_clustering
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

# Figure S12
baker$Total_dissimilarity_plot

# Figure S13
baker$Internal_measures_panel

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

# Figure S14
ggarrange(silh2, silh4, silh8, silh11)


# Table S5: Cophenetic correlation coefficients 
baker$Table_cophenetic_coefficient


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

# Figure S15
ggarrange(distr_char_cluster_inf$`Duration`, 
          distr_char_cluster_inf$`Quality score`, 
          distr_char_cluster_inf$`Inclusion FEV1`, 
          distr_char_cluster_inf$`Inclusion FVC`,
          distr_char_cluster_inf$`Smoking history`, 
          distr_char_cluster_inf$`Minimum FEV1 (%)`, 
          distr_char_cluster_inf$`Maximum FEV1 (%)`, 
          common.legend = TRUE, 
          legend = "bottom")

# Figure S16
ggarrange(distr_char_cluster_inf$`Random allocation`, 
          distr_char_cluster_inf$`Double blinding`, 
          distr_char_cluster_inf$`Withdrawals description`, 
          common.legend = FALSE, 
          legend = "bottom")


## Distribution of characteristics per cluster (HEURISTIC CLUSTERING) ----
distr_char_cluster <- distr_characteristics(input = dataset,
                                            drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                                            rename_char = list(c(5:14),
                                                               c("Duration", "Random allocation", "Double blinding",
                                                                 "Withdrawals description", "Quality score", "Inclusion FEV1",
                                                                 "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                                            cluster = baker,
                                            label_size = 3.0,
                                            axis_text_size = 12)

# Figure S17 
ggarrange(distr_char_cluster$`Duration`, 
          distr_char_cluster$`Quality score`, 
          distr_char_cluster$`Inclusion FEV1`, 
          distr_char_cluster$`Inclusion FVC`,
          distr_char_cluster$`Smoking history`, 
          distr_char_cluster$`Minimum FEV1 (%)`, 
          distr_char_cluster$`Maximum FEV1 (%)`, 
          common.legend = TRUE, 
          legend = "bottom")

# Figure S18
ggarrange(distr_char_cluster$`Random allocation`, 
          distr_char_cluster$`Double blinding`, 
          distr_char_cluster$`Withdrawals description`, 
          common.legend = FALSE, 
          legend = "bottom")
