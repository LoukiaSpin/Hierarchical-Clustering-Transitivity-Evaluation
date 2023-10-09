#*******************************************************************************
#*
#*
#*               Create Main Figures from both datasets
#*               
#*
#*******************************************************************************


## Load rnmamod developmental version
#remotes::install_github("https://github.com/LoukiaSpin/rnmamod.git", force = TRUE)


## Load libraries 
list.of.packages <- c("readxl", "ggplot2", "reshape2", "ggpubr", "heatmaply", "igraph", 
                      "gridExtra", "writexl", "scales", "stringr", "rnmamod")
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


#*******************************************************************************
#*                          SINGH DATASETS & ANALYSES                                                                             
#*******************************************************************************


## Load dataset for clustering - SINGH
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


## Load dataset for network plot - SINGH
#data_network <- as.data.frame(read_excel("./31_Dataset - Extracted networks/2009/2009_19821440_Singh.xlsx", na = "NA"))[, c(1:3, 18:19)]
load("./data/singh_network.RData")

# Remove studies with zero sample size
data_network_singh <- subset(data_network, sample1 > 0)

# Turn treatment columns from 'double' to 'character'
data_network_singh[, 2:3] <- lapply(data_network_singh[, 2:3], as.character)

# Check if typeof has been corrected
lapply(data_network_singh, typeof)

# Turn into the proper format 
network_data_singh <- long_to_wide(input = data_network_singh)

# Conduct *informative* hierarchical clustering
singh_inf <- comp_clustering(input = dataset_singh, 
                             drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                             rename_char = 
                               list(5:13, 
                                    c("Duration", "Disease duration", "MTX use", "RA duration", "Anti-TNF",
                                      "Prior drugs fail", "Prior TNF fail", "Comb. biologic", "Biologic naive")),
                             height = TRUE,
                             get_plots = TRUE,
                             label_size = 4,
                             title_size = 14,
                             axis_title_size = 14,
                             axis_text_size = 14,
                             axis_x_text_angle = 0,
                             legend_text_size = 13)

# Conduct *heuristic* hierarchical clustering
singh_heu <- comp_clustering(input = dataset_singh, 
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



#*******************************************************************************
#*                          BAKER DATASETS & ANALYSES                                                                                                       
#*******************************************************************************


## Load dataset for clustering - BAKER
#dataset_baker <- as.data.frame(read_excel("./31_Dataset - Extracted networks/2009/2009_19637942_Baker.xlsx", na = "NA"))[, -c(4:5, 15, 18:21)]
load("./data/baker_dataset.RData")
dataset_baker <- dataset

# Turn treatment columns from 'double' to 'character'
dataset_baker[, 2:3] <- lapply(dataset_baker[, 2:3], as.character)

# Turn character into factor for the corresponding characteristics
dataset_baker[, 6:8] <- lapply(dataset_baker[, 6:8], as.factor)

# Check if typeof has been corrected
lapply(dataset_baker, typeof)


## Load dataset for network plot - BAKER
#data_network_baker <- as.data.frame(read_excel("./31_Dataset - Extracted networks/2009/2009_19637942_Baker.xlsx", na = "NA"))[, c(1:3, 20:21)]
load("./data/baker_network.RData")
data_network_baker <- data_network

# Turn treatment columns from 'double' to 'character'
data_network_baker[, 2:3] <- lapply(data_network_baker[, 2:3], as.character)

# Check if typeof has been corrected
lapply(data_network_baker, typeof)

# Turn into the proper format 
network_data_baker <- long_to_wide(input = data_network_baker)

# Conduct *informative* hierarchical clustering
baker_inf <- comp_clustering(input = dataset_baker, 
                             drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                             rename_char = list(c(5:14),
                                                c("Duration", "Random allocation", "Double blinding",
                                                  "Withdrawals description", "Quality score", "Inclusion FEV1",
                                                  "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                             height = TRUE,
                             get_plots = TRUE,
                             label_size = 3.0,
                             title_size = 11,
                             axis_title_size = 11,
                             axis_text_size = 9,
                             legend_text_size = 11)

# Conduct *heuristic* hierarchical clustering
baker_heu <- comp_clustering(input = dataset_baker, 
                             drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                             rename_char = list(c(5:14),
                                                c("Duration", "Random allocation", "Double blinding",
                                                  "Withdrawals description", "Quality score", "Inclusion FEV1",
                                                  "Inclusion FVC", "Smoking history", "Minimum FEV1 (%)", "Maximum FEV1 (%)")),
                             height = FALSE,
                             optimal_clusters = 5,
                             get_plots = TRUE,
                             label_size = 3.0,
                             title_size = 11,
                             axis_title_size = 11,
                             axis_text_size = 9,
                             legend_text_size = 11)



#*******************************************************************************
#*                    GET ALL FIGURES (except for Figure 2)                                                                                                                                                                                                                                          
#*******************************************************************************


## Figure 1
tiff("./Figure 1.tiff", 
     height = 20, 
     width = 35, 
     units = 'cm', 
     compression = "lzw", 
     res = 600)
layout(matrix(c(1, 2), ncol =  2, byrow = TRUE))

netplot(data = network_data_singh, 
        drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
        edge_color = singh_inf$Cluster_color,
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
        edge_color = baker_inf$Cluster_color,
        show_multi = TRUE,
        multi_frame = -17,
        node_color = "white", 
        node_frame_color = "grey80",
        node_frame_width = 1.6,
        edge_label_cex = 1,
        direction = FALSE)
mtext(expression(paste(bold("b)"))), side = 3, line = 0,  adj = -0.06)
dev.off()


## Figure 3
tiff("./Figure 3.tiff", 
     height = 18, 
     width = 36, 
     units = 'cm', 
     compression = "lzw", 
     res = 600)
ggarrange(comp_clustering(input = dataset_singh, 
                          drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                          height = FALSE,
                          optimal_clusters = 3,
                          get_plots = TRUE,
                          label_size = 3.0,
                          title_size = 11,
                          axis_title_size = 11,
                          axis_text_size = 9,
                          legend_text_size = 11)$Dissimilarity_comparison,
          comp_clustering(input = dataset_singh, 
                          drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                          height = FALSE,
                          optimal_clusters = 3,
                          get_plots = TRUE,
                          label_size = 3.0,
                          title_size = 11,
                          axis_title_size = 11,
                          axis_text_size = 9,
                          legend_text_size = 11)$Total_dissimilarity_plot,
            labels = c("a)", "b)"),
          common.legend = FALSE,
          legend = "bottom")
dev.off()


## Figure 4
# Save each png/tiff, then bring to pptx and save as tiff
dendro_heatmap(input = singh_heu,
               optimal_dist = "euclidean",
               optimal_link = "average",
               optimal_clusters = 3)

tiff("./Figure 4.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
network_comparisons(data = singh_heu, 
                    node_frame_width = 2,
                    node_label_cex = 1.2,
                    node_label_dist = c(0.0, 2.2, 2.2, 0.0, -2.2, -2.2),
                    edge_level = "all")
legend("bottom", 
       legend = c("Low [0.00, 0.25]", "Moderate (0.25, 0.50]", "High (0.50, 0.75]", "Very high (0.75, 1.00]"), 
       col = c("#A6D854","#E6AB02", "#D95F02", "#E31A1C"), 
       pch = 15,
       text.width = c(0.37, 0.50, 0.39, 0.57), #0.064, 0.083, 0.065, 0.07
       xjust = 0, 
       yjust = 0,
       inset = -0.08,
       xpd = TRUE,
       cex = 1.2,
       bty = "n",
       horiz = TRUE)
dev.off()


## Figure 5
tiff("./Figure 5.tiff", 
     height = 18, 
     width = 36, 
     units = 'cm', 
     compression = "lzw", 
     res = 600)
comp_clustering(input = dataset_baker, 
                drug_names = c("PBO", "BUD", "BUD+", "FLU", "FLU+", "FOR", "SAL", "TIO"),
                height = TRUE,
                optimal_clusters = 5,
                get_plots = TRUE,
                label_size = 3.0,
                title_size = 11,
                axis_title_size = 11,
                axis_text_size = 8.5,
                legend_text_size = 11)$Dissimilarity_comparison
dev.off()


## Figure 6
# Save as png, bring to pptx and save as tiff
dendro_heatmap(input = baker_heu,
               optimal_dist = "euclidean",
               optimal_link = "average",
               optimal_clusters = 5)


## Figure 7
# First, check dimensions
layout(matrix(c(1, 2, 3, 4, 5, 5), ncol = 2, byrow = TRUE), 
       heights = c(0.9, 0.9, 0.1),  
       widths = c(0.9, 0.9))
layout.show(n = 5)
dev.off()

# Then get the plot
tiff("./Figure 7.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
layout(matrix(c(1, 2, 3, 4, 5, 5), ncol =  2, byrow = TRUE), heights = c(1.1, 1.1, 0.2))
par(mar = c(0, 0, 0, 0), oma = c(1.2, 1.2, 1.2, 1.2))

network_comparisons(data = baker_heu, 
                    node_frame_width = 1.5,
                    node_label_cex = 1.5,
                    node_label_dist = 0,
                    edge_level = "low")
mtext(expression(paste(bold("a)"))), side = 3, line = -1.6, adj = 0.20)

network_comparisons(data = baker_heu, 
                    node_frame_width = 1.5,
                    node_label_cex = 1.5,
                    node_label_dist = 0,
                    edge_level = "moderate")
mtext(expression(paste(bold("b)"))), side = 3, line = -1.6,  adj = 0.20)

network_comparisons(data = baker_heu, 
                    optimal_clusters = 5,
                    node_frame_width = 1.5,
                    node_label_cex = 1.5,
                    node_label_dist = 0,
                    edge_level = "high")
mtext(expression(paste(bold("c)"))), side = 3, line = -1.6,  adj = 0.20)

network_comparisons(data = baker_heu, 
                    node_frame_width = 1.5,
                    node_label_cex = 1.5,
                    node_label_dist = 0,
                    edge_level = "very high")
mtext(expression(paste(bold("d)"))), side = 3, line = -1.6,  adj = 0.20)

par(mar = c(0, 0, 0, 0)) # Make the margins 0
plot(1, type = "n", axes = F, xlab = "", ylab = "") 

legend("bottom", 
       legend = c("Low [0.00, 0.25]", "Moderate (0.25, 0.50]", "High (0.50, 0.75]", "Very high (0.75, 1.00]"), 
       col = c("#A6D854","#E6AB02", "#D95F02", "#E31A1C"), 
       pch = 15,
       text.width = c(0.064, 0.083, 0.065, 0.07),
       xjust = 0, 
       yjust = 0,
       inset = 0.2,
       xpd = TRUE,
       cex = 1.7,
       bty = "n",
       horiz = TRUE)
dev.off()
