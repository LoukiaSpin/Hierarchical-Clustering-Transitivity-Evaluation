#*******************************************************************************
#*
#*
#*                      Create Figure 3 of the Manuscript 
#*           (Schematic illustration of the proposed novel approach)                       
#*
#*
#*******************************************************************************


## Load rnmamod developmental version
remotes::install_github("https://github.com/LoukiaSpin/rnmamod.git", force = TRUE)


## Load libraries 
library(rnmamod)


#*******************************************************************************
#*                 CREATE FICTIONAL DATASET FROM SINGH DATASET                                                                                                              
#*******************************************************************************


# Fictional dataset
data_set <- data.frame(Trial_name = as.character(1:7),
                       arm1 = c("1", "1", "1", "1", "1", "2", "2"),
                       arm2 = c("2", "2", "2", "3", "3", "3", "3"),
                       sample = c(140, 145, 150, 40, 45, 75, 80),
                       age = c(18, 18, 18, 48, 48, 35, 35),
                       randomisation = factor(c("yes", "yes", "yes", "no", "no", "no", "no")))

# Check if typeof has been corrected
lapply(data_set, typeof)

# Prepare dataset for the 'gower_distance' function
data_set_gower <- data.frame(trial = as.character(1:dim(data_set)[1]),
                             comp = as.character(paste(data_set$treat2, "vs", data_set$treat1)),
                             data_set[, 4:6])

# Check if typeof has been corrected
lapply(data_set_gower, typeof)

# Apply the 'gower_distance' function
(gower_AB <- gower_distance(input = data_set_gower[1:3, ]))
(gower_AC <- gower_distance(input = data_set_gower[4:5, ]))
(gower_BC <- gower_distance(input = data_set_gower[6:7, ]))

# Apply the 'comp_clustering' function for *Informative* decision
comp_clustering(input = data_set,
                drug_names = c("A", "B", "C"),
                threshold = 0.13,  # General research setting
                informative = TRUE,
                get_plots = TRUE)

# Apply the 'comp_clustering' function for *Heuristic* clustering
data_clust <- comp_clustering(input = data_set,
                              drug_names = c("A", "B", "C"),
                              informative = FALSE,
                              optimal_clusters = 3,
                              get_plots = TRUE)

# Get dendrogram
dendro_heatmap(input = data_clust, 
               label_size = 14, 
               axis_text_size = 14)
