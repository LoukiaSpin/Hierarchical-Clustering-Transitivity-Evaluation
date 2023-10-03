#' Dendrogram with amalgamated heatmap
#' (comparisons' comparability for transitivity evaluation)
#' 
#' @description 
#'   \code{dendro_heatmap} creates a dendrogram alongside the heatmap of
#'   dissimilarities among the comparisons for a specific linkage method and
#'   number of clusters.
#'   
#' @param input An object of S3 class \code{\link{comp_clustering}}. See 'Value'
#'   in \code{\link{comp_clustering}}.
#' @param optimal_dist A character string with values \code{"euclidean"}, 
#'   \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, 
#'   \code{"minkowski"}, or \code{"gower"} for the optimal
#'   dissimilarity measure, corresponding to the highest cophenetic correlation
#'   coefficient. See 'Details'.
#' @param optimal_link A character string with values \code{"ward.D"}, 
#'   \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"},
#'   \code{"mcquitty"}, \code{"median"}, or \code{"centroid"} for the optimal
#'   linkage method, corresponding to the highest cophenetic correlation
#'   coefficient. See 'Details'.
#' @param optimal_clusters A positive integer for the optimal number of clusters 
#'   based on three internal validation measures. It takes values from two to 
#'   the number of comparisons minus one. See 'Details'.
#'
#' @return 
#'   \code{dendro_heatmap} uses the \code{\link[heatmaply:heatmaply]{heatmaply}} 
#'   function of the R-package 
#'   \href{https://CRAN.R-project.org/package=heatmaply}{heatmaply} to create a 
#'   cluster heatmap for a selected linkage method and number of clusters. The
#'   function uses different colours to indicate the clusters directly on the
#'   dendrogram.
#'  
#'  @details 
#'    The user needs first to inspect the results of three internal measures 
#'    (connectivity index, silhouette width, and Dunn index) for a wide range of 
#'    clusters to define the argument \code{optimal_clusters}. The user also
#'    needs to inspect the cophenetic correlation coefficient for all pairwise 
#'    combinations of six dissimilarity measures (Euclidean, maximum, Manhattan, 
#'    Canberra, Minkowski, and Gower) with eight linkage methods (two Ward 
#'    versions, single, complete, average, Mcquitty, median and centroid). All 
#'    these inspections are performed using the 
#'    \code{\link{internal_measures_plot}} and \code{\link{comp_clustering}} 
#'    functions, respectively.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso 
#'  \code{\link{comp_clustering}}, \code{\link[heatmaply:heatmaply]{heatmaply}},
#'  \code{\link{internal_measures_plot}}
#'
#' @export
dendro_heatmap <- function (input, 
                            optimal_dist,
                            optimal_link, 
                            optimal_clusters) {
  
  
  ## Check the defaults
  # Dataset
  diss <- if (class(input$Dissimilarity_table) != "dist") {
    stop("'input' must be of class 'dist'", call. = FALSE)
  } else {
    input$Dissimilarity_table
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
  link_list <- c("ward.D", "ward.D2", "single", "complete", "average", 
                    "mcquitty", "median", "centroid")
  m_list1 <- c("'ward.D', 'ward.D2', 'single', 'complete', 'average'")
  m_list2 <- c("'mcquitty', 'median', 'centroid'")
  optimal_link <- if (missing(optimal_link)) {
    stop("The argument 'optimal_link' must be defined", call. = FALSE)
  } else if (!is.element(optimal_link, link_list)) {
    stop(paste("'optimal_link' must be any of the following:", 
               m_list1, m_list2), call. = FALSE)
  } else {
    optimal_link
  }

  # Number of 'optimal' clusters (based on the internal measures)
  compar <- colnames(as.matrix(diss))
  optimal_clusters <- if (missing(optimal_clusters)) {
    2
    #stop("The argument 'optimal_clusters' must be defined", call. = FALSE)
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

  
  ## Function for first letter capital (Source: https://stackoverflow.com/questions/18509527/first-letter-to-upper-case)
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  # Turn first letter capital
  optimal_dist_new <- if (is.element(optimal_dist, 
                                     c(dist_list[-2], "gower"))) {
    firstup(optimal_dist)
  } else {
    optimal_dist
  }

  
  ## Create heatmap with dendrogram and coloured clusters
  #library(heatmaply)
  dendro_heatmap <- 
    heatmaply(as.matrix(diss),
              cellnote = round(as.matrix(diss), 2),
              xlab = " ",
              ylab = " ",
              scale = "none",
              Colv = reorder(as.dendrogram(hclust(diss)), 
                             input$Total_dissimilarity[, 3], max), 
              Rowv = reorder(as.dendrogram(hclust(diss)), 
                             input$Total_dissimilarity[, 3], max), 
              k_col = optimal_clusters,
              k_row = optimal_clusters,
              scale_fill_gradient_fun = 
                scale_fill_gradient2(name = paste(optimal_dist_new, 
                                                  "dissimilarity"), 
                                     low = "white", 
                                     high = "red", 
                                     na.value = "grey90"#,
                                     #limit = c(limits_scale[1], 
                                     #          limits_scale[2]))
              ))  
  
  
  return(dendro_heatmap)
} 