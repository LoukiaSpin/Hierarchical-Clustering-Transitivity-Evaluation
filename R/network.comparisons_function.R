#' Network of pairwise comparisons of interventions
#' (Comparisons' comparability for transitivity evaluation)
#' 
#' @description 
#'   \code{network_comparisons} offers an  alternative visualisation of the 
#'   dendrogram with amalgamated heatmap. Presenting the comparisons in a 
#'   circular layout, typically seen in networks, facilitates interpretation of
#'   the clustering results in the context of transitivity evaluation.
#'   
#' @param data A data-frame in the long arm-based format. Two-arm trials occupy
#'   one row in the data-frame. Multi-arm trials occupy as many rows as the
#'   number of possible comparisons among the interventions. The first three 
#'   columns refer to the trial name, first and second arm of the comparison, 
#'   respectively. The remaining columns refer to summary characteristics. See 
#'   'Details' for the specification of the columns. This is the same format
#'   required to apply \code{\link{comp_clustering}}, the core function to 
#'   conduct clustering of comparisons.
#' @param drug_names A vector of labels with the name of the interventions
#'   in the order they appear in the argument \code{input}.
#' @param total_diss An object of S3 class \code{\link{comp_clustering}} with
#'   the total dissimilarities of the observed comparisons. See 'Value' in 
#'   \code{\link{comp_clustering}}.
#' @param comp_diss An object of S3 class \code{\link{comp_clustering}} with
#'   the dissimilarities among the observed comparisons. See 'Value' in 
#'   \code{\link{comp_clustering}}.
#' @param optimal_link A character string with values \code{"ward.D"}, 
#'   \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"},
#'   \code{"mcquitty"}, \code{"median"}, or \code{"centroid"} for the optimal
#'   linkage method, corresponding to the highest cophenetic correlation
#'   coefficient. See 'Details' below.
#' @param optimal_clusters A positive integer for the optimal number of clusters 
#'   based on three internal validation measures. It takes values from two to 
#'   the number of comparisons minus one. 
#' @param node_frame_color A character string for the colour of the node frames. 
#'   The default argument is "black". For more details refer to the R-package 
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}. 
#' @param node_frame_width A positive integer for the size of the node frame.
#'   The default argument is 1. For more details refer to the R-package 
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param node_label_font A positive integer for the font type of the node 
#'   labels. The default argument is 1. For more details refer to the R-package 
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param node_label_cex A positive integer for the size of the node labels.
#'   The default argument is 1. For more details refer to the R-package 
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param node_label_dist A positive integer for the distance between label and
#'   node. The default argument is 0. For more details refer to the R-package 
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param node_label_color A character string for the colour of the node labels. 
#'   The default argument is "black". For more details refer to the R-package 
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param edge_level A character string with values \code{"low"}, \code{"low"}, 
#'   \code{"moderate"}, \code{"high"}, and \code{"very high"} to indicate the 
#'   level of dissimilarity between two comparisons to draw. 
#'   See 'Details' below.
#' @param edge_lty A positive integer for the line type of the edges. The 
#'   default argument is 1. For more details refer to the R-package 
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param edge_curved A value in the range 0 to 1 for the edge curvature. The 
#'   default argument is 0. For more details refer to the R-package 
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param ... Further graphical arguments of the 
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}
#'   
#' @return
#'   The network of comparisons with different colours for the clustered nodes.
#'
#' @details
#'    The nodes refer to the observed comparisons and their size is proportional
#'    to the corresponding normalised total dissimilarity. The nodes are 
#'    coloured with different colours to indicate the cluster they belong. 
#'    The clusters are determined by the cutree function for specific 
#'    \code{optimal_link} and \code{optimal_clusters}. The edges refer to the 
#'    dissimilarities of pairs of comparisons and their size is proportional to 
#'    the corresponding normalised dissimilarity value. Four different colours 
#'    are used to indicate the level of dissimilarity between two comparisons:
#'    green, yellow, orange, and red for low, moderate, high and very high 
#'    dissimilarity (\eqn{d < 0.25}, \eqn{0.25 \qeq d < 0.5}, 
#'    \eqn{0.50 \qeq d < 0.75}, and \eqn{d \qeq 0.75}, respectively).
#'    
#'    The interventions should be sorted in an ascending order of their 
#'    identifier number within the trials so that the first treatment column 
#'    (second column in \code{data}) is the control arm for every pairwise 
#'    comparison. This is important to ensure consistency in the order of 
#'    interventions within the comparisons obtained from the other related 
#'    functions and to colour the correct comparisons via the 
#'    \code{\link{netplot}} function.
#'
#'    Networks with many comparisons appear cluttered when showing all edges
#'    (i.e., \code{edge_level = "all"}). In this case, it is recommended to draw 
#'    the network for the dissimilarity level of interest, for instance,
#'    \code{edge_level = "low"}). When a dissimilarity level does not exist, the 
#'    network will be drawn without the corresponding edges.
#'    
#'    The user needs first to inspect the results of three internal measures 
#'    (connectivity index, silhouette width, and Dunn index) for a wide range of 
#'    clusters to define the argument \code{optimal_clusters}. The user also
#'    needs to inspect the cophenetic correlation coefficient for all pairwise 
#'    combinations of six dissimilarity measures (Euclidean, maximum, Manhattan, 
#'    Canberra, Minkowski, and Gower) with eight linkage methods (two Ward 
#'    versions, single, complete, average, McQuitty, median and centroid). All 
#'    these inspections are performed using the 
#'    \code{\link{internal_measures_plot}} and \code{\link{comp_clustering}} 
#'    functions, respectively.
#' 
#' @author {Loukia M. Spineli}
#'
#' @seealso 
#'  \code{\link{comp_clustering}}, \code{\link{internal_measures_plot}}, 
#'  \code{\link{netplot}}, \code{\link[igraph:plot.igraph]{plot.igraph}}
#'
#' @export
network_comparisons <- function(data, 
                                drug_names,
                                total_diss, 
                                comp_diss,
                                optimal_dist,
                                optimal_link,
                                optimal_clusters,
                                node_frame_color = "black",
                                node_frame_width = 1,
                                node_label_font = 1,
                                node_label_cex = 1,
                                node_label_dist = 0, 
                                node_label_color = "black",
                                edge_level,
                                edge_lty = 1,
                                #edge_label_color = "black",
                                #edge_label_font = 1,
                                #edge_label_cex = 1,
                                edge_curved = 0,
                                ...) {
  
  
  ## Check the defaults
  # Dataset
  data <- if (any(sapply(data, typeof)[1:3] != "character")) {
    stop("The first three columns (trial and arms) must be 'characters'.", 
         call. = FALSE)
  } else if (any(sapply(data, typeof)[-c(1:3)] == "character")) {
    stop("The characteristics must be 'double' or 'integer'.", call. = FALSE)
  } else {
    data
  }
  colnames(data)[1:3] <- c("Trial_name", "Arm1", "Arm2")
  
  # Intervention names
  drug_names <- if (missing(drug_names)) {
    1:length(unique(data[, 2:3]))
  } else {
    drug_names
  }
  
  # 'Optimal' dissimilarity method (based on the cophenetic coefficient)
  dist_list <- 
    c("euclidean", "maximum", "manhattan", "canberra", "minkowski", "gower")
  d_list <- 
    c("'euclidean', 'maximum', 'manhattan', 'canberra', 'minkowski', 'gower'")
  optimal_dist <- if (missing(optimal_dist)) {
    stop("The argument 'optimal_dist' must be defined", call. = FALSE)
  } else if (!is.element(optimal_dist, dist_list)) {
    stop(paste("'optimal_dist' must be any of the following:", d_list), 
         call. = FALSE)
  } else {
    optimal_dist
  }
  
  # 'Optimal' linkage method (based on the cophenetic coefficient)
  methods_list <- c("ward.D", "ward.D2", "single", "complete", "average", 
                    "mcquitty", "median", "centroid")
  m_list1 <- c("'ward.D', 'ward.D2', 'single', 'complete', 'average'")
  m_list2 <- c("'mcquitty', 'median', 'centroid'")
  optimal_link <- if (missing(optimal_link)) {
    stop("The argument 'optimal_link' must be defined", call. = FALSE)
  } else if (!is.element(optimal_link, methods_list)) {
    stop(paste("'optimal_link' must be any of the following:", 
               m_list1, m_list2), call. = FALSE)
  } else {
    optimal_link
  }
  
  # Number of 'optimal' clusters (based on the internal measures)
  compar <- colnames(as.matrix(data))
  optimal_clusters <- if (missing(optimal_clusters)) {
    2
  } else if ((optimal_clusters > length(compar) - 1 ||
              optimal_clusters < 2) & length(compar) > 3) {
    stop(paste0("'optimal_clusters' must range from 2 to", " ", 
                length(compar) - 1, "."), call. = FALSE)
  } else if ((optimal_clusters > length(compar) - 1 ||
              optimal_clusters < 2) & length(compar) == 3) {
    stop(paste0("'optimal_clusters' must equal exactly 2."), call. = FALSE)
  } else {
    optimal_clusters
  }
  
  # Set the default for 'edge_level'
  levels <- c("low", "moderate", "high", "very high", "all")
  edge_level <- if (missing(edge_level)) {
    a <- c("'low', 'moderate', 'high', 'very high', 'all'")
    stop(paste("The argument 'edge_level' must be one of the following:", a))
  } else if (!is.element(edge_level, levels)) {
    a <- c("'low', 'moderate', 'high', 'very high', 'all'")
    stop(paste("The argument 'edge_level' must be one of the following:", a))
  } else {
    edge_level
  }
  
  
  ## Assign the intervention names (if applicable)
  data[, 2:3] <- matrix(drug_names[as.numeric(unlist(data[, 2:3]))], 
                        nrow = dim(data)[1],
                        ncol = 2)

  
  ## Extract the columns with the intervention id
  pairwise <- data[, startsWith(colnames(data), "Arm")]
  
  
  ## Name the comparisons (control appears second in each comparison)
  pairwise_name <- unique(paste(pairwise[, 2], "vs", pairwise[, 1]))

  
  ## Possible pairwise comparisons (reflect matrix of total dissimilarities)
  poss_comp <- combn(pairwise_name, 2)

  
  ## Prepare plot (igraph package)
  g1 <- igraph::graph(edges = poss_comp, directed = FALSE) 

  
  ## Sort comparisons (nodes) by the order of 'poss_comp'
  total_diss_new <- 
    total_diss[match(unique(c(poss_comp)), total_diss$comparison), ]
  
  
  ## Weight each node by the corresponding total distance (normalised)
  #V(g1)$weight <- 
  #  (0.40 + ((total_diss_new[, 3] - min(total_diss_new[, 3])) / 
  #             (max(total_diss_new[, 3]) - min(total_diss_new[, 3])) )) * 20
  
  
  ## Weight each node by the corresponding total distance
  igraph::V(g1)$weight <- (0.40 + total_diss_new[, 3]) * 20
  
  
  ## Match comparisons with their cluster
  cut_dend <- cutree(hclust(comp_diss, method = optimal_link),
                     k = optimal_clusters) 

  
  ## Turn into data-frame
  cluster_comp0 <- data.frame(comp = rownames(data.frame(cut_dend)),
                              cluster = cut_dend)
  rownames(cluster_comp0) <-  NULL
  
  
  ## Sort comparisons by the order of 'poss_comp'
  cluster_comp <- 
    cluster_comp0[match(unique(c(poss_comp)), cluster_comp0$comp), ]
  
  
  ## Colour node by cluster
  #node_col <- rev(hue_pal()(optimal_clusters))[cluster_comp[, 2]] 
  node_col <- hue_pal()(optimal_clusters)[cluster_comp[, 2]] 

  
  ## Possible comparisons among comparisons
  poss_comp_comp <- combn(pairwise_name, 2)
  
  
  ## Sort comparisons of comparisons (edges) by the order of 'poss_comp'
  comp_diss_new <- rep(NA, dim(poss_comp_comp)[2])
  for (i in 1:dim(poss_comp_comp)[2]) {
    comp_diss_new[i] <- as.matrix(comp_diss)[poss_comp_comp[1,i], poss_comp_comp[2,i]]
  }

  
  ## Normalise the dissimilarity matrix of comparisons
  norm_comp_diss <- if (optimal_dist == "gower") {
    comp_diss_new
  } else {
    (comp_diss_new - min(comp_diss_new)) / 
      (max(comp_diss_new) - min(comp_diss_new))
  }


  ## Weight each edge by the corresponding distance (normalised)
  igraph::E(g1)$weight <- (0.40 + norm_comp_diss) * 10
  
  
  ## Name the edges according to dissimilarity value
  edge_col_num <- 
    as.vector(do.call(rbind, 
                      lapply(norm_comp_diss, 
                             function(x) if (x < 0.25) "low" 
                             else if (x >= 0.25 & x < 0.50) "moderate" 
                             else if (x >= 0.50 & x < 0.75) "high"
                             else "very high")))
  
  
  ## Assign a colour to the name
  edge_colour <- 
    as.vector(do.call(rbind, 
                      lapply(edge_col_num, 
                             function(x) if (x == "low") "#A6D854" 
                             else if (x == "moderate") "#E6AB02"
                             else if (x == "high") "#D95F02"
                             else "#E31A1C")))
  
  
  ## Adjust the selected level for larger networks (at least 5 comparisons)
  edge_col_select0 <- subset(edge_colour, edge_col_num == edge_level)
  edge_col_select <- ifelse(!is.element(edge_colour, edge_col_select0), 
                            NA, 
                            edge_col_select0)
  
  
  ## Colour the edges based on the network size
  if (edge_level != "all") {
    
    # Adjust the selected level for larger networks (at least 5 comparisons)
    edge_col0 <- subset(edge_colour, edge_col_num == edge_level)
    edge_col <- ifelse(!is.element(edge_colour, edge_col0), NA, edge_col0)
  } else {
    edge_col <- edge_colour
  }

  
  ## Colour edge by cluster
  # Indicate the edges that belong to the same cluster
  #same_cluster <- mapply(function(x, y) setequal(x, y),
  #                       t(combn(cut_dend, 2))[, 1],
  #                       t(combn(cut_dend, 2))[, 2])
  
  # Replace the edge with the cluster number
  #edge_col_num <- 
  #  ifelse(same_cluster == TRUE, apply(t(combn(cut_dend, 2)), 1, min),
  #         apply(t(combn(cut_dend, 2)), 1, max))
  
  # Assign a colour to the number
  #edge_col <- colours[edge_col_num]
  #edge_col <- c("grey80", "black", "orange")[edge_col_num]

  
  ## Label the edges 
  #E(g1)$names <- round(as.vector(comp_diss), 2)
  
  
  ## Create the igraph
  plot(g1,
       layout = layout_in_circle,
       vertex.color = node_col,   
       vertex.frame.color = node_frame_color,
       vertex.frame.width = node_frame_width,
       vertex.shape = "circle",
       vertex.size = igraph::V(g1)$weight,
       vertex.label.font = node_label_font,
       vertex.label.cex	= node_label_cex,
       vertex.label.dist = node_label_dist,
       vertex.label.color = node_label_color,
       edge.color = edge_col, 
       edge.width = igraph::E(g1)$weight,
       edge.lty = edge_lty,
       #edge.label = edge_name, #E(g1)$names, 
       #edge.label.color = edge_label_color,
       #edge.label.font = edge_label_font,
       #edge.label.cex = edge_label_cex,
       edge.curved = edge_curved,
       ...)
  
  
  ## Add legend
  #legend('bottomleft',
  #       title = "Comparison dissimilarity (normalised)",
  #       legend = c("Low [0.00, 0.25)", 
  #                  "Moderate [0.25, 0.50)", 
  #                  "High [0.50, 0.75)", 
  #                  "Very high [0.75, 1.00]"),
  #       fill = c("#A6D854", "#E6AB02", "#D95F02", "#E31A1C"),
  #       inset = .02,
  #       cex = 1.2,
  #       bty = "n",
  #       horiz = FALSE)
  
}