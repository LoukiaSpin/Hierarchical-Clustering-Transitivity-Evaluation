#*******************************************************************************
#*
#*  R Code to replicate the Supplementary Figures & Tables 
#*  <Singh and colleagues (2009; PMID: 19821440) - Rheumatoid arthritis network>
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


## Load Singh dataset for clustering ----
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
# ?distr_characteristics
distr_char <- distr_characteristics(input = dataset_singh,
                                    drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                                    rename_char = 
                                      list(5:13, 
                                           c("Duration", "Disease duration", "MTX use", "RA duration", "Anti-TNF",
                                             "Prior drugs fail", "Prior TNF fail", "Comb. biologic", "Biologic naive")),
                                    label_size = 5,
                                    axis_text_size = 11,
                                    axis_x_text_angle = 0)

# Figure S1
ggarrange(distr_char$Duration, 
          distr_char$`Disease duration`,
          common.legend = TRUE, 
          legend = "bottom")

# Figure S2
ggarrange(distr_char$`MTX use`, 
          distr_char$`Anti-TNF`, 
          distr_char$`Prior drugs fail`,
          distr_char$`Comb. biologic`, 
          distr_char$`Biologic naive`,
          distr_char$`RA duration`, 
          distr_char$`Prior TNF fail`, 
          common.legend = FALSE, 
          legend = "bottom")


## Conduct *HEURISTIC* clustering of comparisons ----
# ?comp_clustering
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

# Figure S3
singh$Internal_measures_panel

# Figure S4
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

# Table S4: Cophenetic correlation coefficients 
singh$Table_cophenetic_coefficient


## Conduct *INFORMATIVE* clustering of comparisons ----
singh_inf <- comp_clustering(input = dataset_singh, 
                             drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),  
                             rename_char = 
                               list(5:13, 
                                    c("Duration", "Disease duration", "MTX use", "RA duration", "Anti-TNF",
                                      "Prior drugs fail", "Prior TNF fail", "Comb. biologic", "Biologic naive")),
                             height = TRUE,
                             get_plots = TRUE)


## Distribution of characteristics per cluster ----
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

# Figure S5
ggarrange(distr_char_cluster$Duration, 
          distr_char_cluster$`Disease duration`,
          common.legend = TRUE, 
          legend = "bottom")

# Figure 6
ggarrange(distr_char_cluster$`MTX use`, 
          distr_char_cluster$`Anti-TNF`, 
          distr_char_cluster$`Prior drugs fail`,
          distr_char_cluster$`Comb. biologic`, 
          distr_char_cluster$`Biologic naive`,
          distr_char_cluster$`RA duration`, 
          distr_char_cluster$`Prior TNF fail`,
          common.legend = FALSE, 
          legend = "bottom")
