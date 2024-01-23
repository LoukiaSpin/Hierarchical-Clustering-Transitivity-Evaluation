#*******************************************************************************
#*
#*  R Code to replicate the Supplementary Figures & Tables 
#*  <Singh and colleagues (2009; PMID: 19821440) - Rheumatoid arthritis network>
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
                                    label_size = 3,
                                    title_size = 9.5,
                                    axis_title_size = 9.5,
                                    axis_text_size = 9.5,
                                    axis_x_text_angle = 0,
                                    legend_text_size = 9.5)

# Figure S1
tiff("./Figure S1.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
ggarrange(distr_char$Duration, 
          distr_char$`Disease duration`,
          common.legend = TRUE, 
          legend = "bottom")
dev.off()

# Figure S2
tiff("./Figure S2.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
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


## Create dissimilarity matrix (D) and heatmap ----
# ?comp_clustering
singh <- comp_clustering(input = dataset_singh, 
                         drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),  
                         threshold = 0.28, #mean(dataset_singh$sample.size) -> 283.0741
                         informative = TRUE,
                         get_plots = TRUE,
                         label_size = 4,
                         title_size = 14,
                         axis_title_size = 14,
                         axis_text_size = 12,
                         axis_x_text_angle = 0,
                         legend_text_size = 13,
                         str_wrap_width = 5)


## Conduct hierarchical clustering ----
# ?comp_clustering
singh_hie <- comp_clustering(input = dataset_singh, 
                             drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),  
                             informative = FALSE,
                             optimal_clusters = 2,
                             get_plots = TRUE,
                             label_size = 4,
                             title_size = 14,
                             axis_title_size = 14,
                             axis_text_size = 12,
                             axis_x_text_angle = 0,
                             legend_text_size = 13,
                             str_wrap_width = 5)

# Table S4: Cophenetic correlation coefficients 
singh_hie$Table_cophenetic_coefficient

# Figure S3
tiff("./Figure S3.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
singh$Within_comparison_dissimilarity
dev.off()

# Figure S4
tiff("./Figure S4.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
singh$Between_comparison_dissimilarity
dev.off()

# Figure S5
tiff("./Figure S5.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
singh_hie$Profile_plot
dev.off()

# Figure S6
tiff("./Figure S6.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
ggarrange(comp_clustering(input = dataset_singh, 
                          drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),  
                          informative = FALSE,
                          optimal_clusters = 2,
                          get_plots = TRUE,
                          label_size = 2.5,
                          title_size = 10,
                          axis_title_size = 10,
                          axis_text_size = 9,
                          legend_text_size = 10)$Silhouette_width_plot,
          comp_clustering(input = dataset_singh, 
                          drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),  
                          informative = FALSE,
                          optimal_clusters = 3,
                          get_plots = TRUE,
                          label_size = 2.5,
                          title_size = 10,
                          axis_title_size = 10,
                          axis_text_size = 9,
                          legend_text_size = 10)$Silhouette_width_plot,
          comp_clustering(input = dataset_singh, 
                          drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),  
                          informative = FALSE,
                          optimal_clusters = 4,
                          get_plots = TRUE,
                          label_size = 2.5,
                          title_size = 10,
                          axis_title_size = 10,
                          axis_text_size = 9,
                          legend_text_size = 10)$Silhouette_width_plot,
          comp_clustering(input = dataset_singh, 
                          drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),  
                          informative = FALSE,
                          optimal_clusters = 5,
                          get_plots = TRUE,
                          label_size = 2.5,
                          title_size = 10,
                          axis_title_size = 10,
                          axis_text_size = 9,
                          legend_text_size = 10)$Silhouette_width_plot,
          ncol = 2,
          nrow = 2,
          labels = c("a)", "b)", "c)", "d)"),
          font.label = list(size = 11))
dev.off()

# Figure S7
dendro_heatmap(input = singh_hie) 


## Distribution of characteristics per cluster ----
distr_char_cluster <- distr_characteristics(input = dataset_singh,
                                            drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                                            rename_char = 
                                              list(5:13, 
                                                   c("Duration", "Disease duration", "MTX use", "RA duration", "Anti-TNF",
                                                     "Prior drugs fail", "Prior TNF fail", "Comb. biologic", "Biologic naive")),
                                            cluster = singh_hie,
                                            label_size = 2.5,
                                            title_size = 9.5,
                                            axis_title_size = 9.5,
                                            axis_text_size = 9.5,
                                            axis_x_text_angle = 0,
                                            legend_text_size = 9.5)

# Figure S8
tiff("./Figure S8.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
ggarrange(distr_char_cluster$Duration, 
          distr_char_cluster$`Disease duration`,
          common.legend = TRUE, 
          legend = "bottom")
dev.off()

# Figure S9
tiff("./Figure S9.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
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