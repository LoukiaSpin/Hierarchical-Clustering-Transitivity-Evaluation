#*******************************************************************************
#*
#*  R Code to replicate the Main Figures
#*  <Singh and colleagues (2009; PMID: 19821440) - Rheumatoid arthritis network>
#*  <Baker and colleagues (2009; PMID: 19637942) - COPD network>
#*     
#*  Authors: Loukia M. Spineli, Katerina Papadimitropoulou, Chrysostomos Kalyvas
#*  Date: December 2023
#*  
#*******************************************************************************


## Load rnmamod developmental version
remotes::install_github("https://github.com/LoukiaSpin/rnmamod.git", force = TRUE)


## Load libraries ----
list.of.packages <- c("rnmamod", "reshape2", "ggpubr", "stringr", "cluster", "scales")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)


## Load functions ----
source("./R/long.to.wide_function.R")


#*******************************************************************************
#*
#*           Data & analysis for Singh (Rheumatoid arthritis network)                                                                                                          
#* 
#*******************************************************************************


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


## Load Singh dataset for network plot ----
load("./data/singh_network.RData")

# Remove studies with zero sample size
data_network_singh <- subset(data_network, sample1 > 0)

# Turn treatment columns from 'double' to 'character'
data_network_singh[, 2:3] <- lapply(data_network_singh[, 2:3], as.character)

# Check if typeof has been corrected
lapply(data_network_singh, typeof)

# Turn into the proper format 
network_data_singh <- long_to_wide(input = data_network_singh)


## Create dissimilarity matrix (D) and heatmap (Singh) ----
# ?comp_clustering
singh_inf <- comp_clustering(input = dataset_singh, 
                             drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                             threshold = 0.28, #mean(dataset_singh$sample.size) -> 283.0741
                             informative = TRUE,
                             get_plots = TRUE,
                             label_size = 4,
                             title_size = 14,
                             axis_title_size = 14,
                             axis_text_size = 14,
                             axis_x_text_angle = 0,
                             legend_text_size = 13,
                             str_wrap_width = 5)


## Conduct hierarchical clustering (Singh) ----
# ?comp_clustering
singh_hie <- comp_clustering(input = dataset_singh, 
                             drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                             informative = FALSE,
                             optimal_clusters = 2,
                             get_plots = TRUE,
                             label_size = 4,
                             title_size = 14,
                             axis_title_size = 14,
                             axis_text_size = 14,
                             axis_x_text_angle = 0,
                             legend_text_size = 13,
                             str_wrap_width = 8)


#*******************************************************************************
#*
#*                  Data & analysis for Baker (COPD network)                                                                                                                                           
#* 
#*******************************************************************************


## Load Baker dataset for clustering ----
load("./data/baker_dataset.RData")
dataset_baker <- dataset

# Turn treatment columns from 'double' to 'character'
dataset_baker[, 2:3] <- lapply(dataset_baker[, 2:3], as.character)

# Turn character into factor for the corresponding characteristics
dataset_baker[, 6:8] <- lapply(dataset_baker[, 6:8], as.factor)

# Check if typeof has been corrected
lapply(dataset_baker, typeof)


## Load Baker dataset for network plot ----
load("./data/baker_network.RData")
data_network_baker <- data_network

# Turn treatment columns from 'double' to 'character'
data_network_baker[, 2:3] <- lapply(data_network_baker[, 2:3], as.character)

# Check if typeof has been corrected
lapply(data_network_baker, typeof)

# Turn into the proper format 
network_data_baker <- long_to_wide(input = data_network_baker)


## Create dissimilarity matrix (D) and heatmap (Baker) ----
# ?comp_clustering
baker_inf <- comp_clustering(input = dataset_baker, 
                             weight = c(rep(1, 9), 0.5, 0.5),
                             drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                             threshold = 0.13,
                             informative = TRUE,
                             get_plots = TRUE,
                             label_size = 3.0,
                             title_size = 11,
                             axis_title_size = 11,
                             axis_text_size = 9,
                             legend_text_size = 11)


## Conduct hierarchical clustering (Baker) ----
# ?comp_clustering
baker_hie <- comp_clustering(input = dataset_baker, 
                             weight = c(rep(1, 9), 0.5, 0.5),
                             drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                             informative = FALSE,
                             optimal_clusters = 2, 
                             get_plots = TRUE,
                             label_size = 3.0,
                             title_size = 11,
                             axis_title_size = 11,
                             axis_text_size = 9.8,
                             legend_text_size = 11)


## Figure 1
tiff("./Figure 1.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
layout(matrix(c(1, 2), ncol =  2, byrow = TRUE))

netplot(data = network_data_singh, 
        drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
        show_multi = TRUE,
        multi_frame = -14,
        node_color = "white", 
        node_frame_color = "grey80",
        node_frame_width = 1.6,
        edge_label_cex = 1,
        direction = FALSE)
mtext(expression(paste(bold("a)"))), side = 3, line = 0,  adj = -0.06)

netplot(data = network_data_baker, 
        drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
        show_multi = TRUE,
        multi_frame = -17,
        node_color = "white", 
        node_frame_color = "grey80",
        node_frame_width = 1.6,
        edge_label_cex = 1,
        direction = FALSE)
mtext(expression(paste(bold("b)"))), side = 3, line = 0,  adj = -0.06)
dev.off()


## Figure 4
tiff("./Figure 4.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
singh_inf$Dissimilarity_heatmap
dev.off()


## Figure 5
tiff("./Figure 5.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
singh_hie$Barplot_comparisons_cluster
dev.off()


## Figure 6
tiff("./Figure 6.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
baker_inf$Dissimilarity_heatmap
dev.off()


## Figure 7
tiff("./Figure 7.tiff", height = 20, width = 35, units = 'cm', compression = "lzw", res = 600)
baker_hie$Barplot_comparisons_cluster
dev.off()

