#' Perform clustering within a level set
#'
#' @param points_in_this_level Points in the current level set.
#' @param filter_values The filter values.
#' @param num_bins_when_clustering Number of bins when clustering.
#' @param methods Specify the clustering method to be used, e.g., "hclust" or "kmeans".
#' @param max_kmeans_clusters Maximum number of clusters when using k-means clustering.
#' @param eps The maximum distance between two samples for one to be considered as in the neighborhood of the other.
#' @param minPts The number of samples in a neighborhood for a point to be considered as a core point.
#' @param num_clusters Number of clusters when using PAM clustering.
#' @return A list containing the number of vertices, external indices, and internal indices.
#' @importFrom stats as.dist hclust cutree dist kmeans
#' @export
perform_clustering <- function(
    points_in_this_level, filter_values, num_bins_when_clustering, methods, 
    max_kmeans_clusters = 10,  # Kmeans
    eps = 0.5, minPts = 5, # DBSCAN
    num_clusters = 5 # PAM
    ) {
  num_points_in_this_level <- length(points_in_this_level)
  
  if (num_points_in_this_level == 0) {
    return(list(num_vertices = 0, external_indices = NULL, internal_indices = NULL))
  }
  
  if (num_points_in_this_level == 1) {
    return(list(num_vertices = 1, external_indices = points_in_this_level, internal_indices = c(1)))
  }

  level_dist_object <- as.dist(as.matrix(dist(filter_values))[points_in_this_level, points_in_this_level])

  if (methods == "hierarchical"){
    level_max_dist <- max(level_dist_object)
    level_hclust <- hclust(level_dist_object, method = "single")
    level_heights <- level_hclust$height

    level_cutoff <- cluster_cutoff_at_first_empty_bin(level_heights, level_max_dist, num_bins_when_clustering)
    level_external_indices <- points_in_this_level[level_hclust$order]
    level_internal_indices <- as.vector(cutree(list(
      merge = level_hclust$merge,
      height = level_hclust$height,
      labels = level_external_indices), h = level_cutoff))
    num_vertices_in_this_level <- max(level_internal_indices)
    
  } else if (methods == "kmeans") {
    max_clusters_for_this_level <- min(max_kmeans_clusters, num_points_in_this_level)
    
    level_filter_values <- filter_values[points_in_this_level, , drop = FALSE]
  
    # If there are more than one point in the level, perform k-means clustering
    if (max_clusters_for_this_level < nrow(level_filter_values)) {
      level_kmean <- kmeans(level_filter_values, centers = max_clusters_for_this_level)
      
      level_external_indices <- points_in_this_level[order(level_kmean$cluster)]
      level_internal_indices <- as.vector(level_kmean$cluster)
      num_vertices_in_this_level <- max(level_internal_indices)
    } else {
      level_external_indices <- points_in_this_level
      level_internal_indices <- rep(1, num_points_in_this_level)
      num_vertices_in_this_level <- 1
    }
  } else if (methods == "dbscan") {
    level_filter_values <- filter_values[points_in_this_level, , drop = FALSE]

    dbscan_result <- dbscan::dbscan(level_filter_values, eps = eps, minPts = minPts)

    # If DBSCAN finds clusters
    if (max(dbscan_result$cluster) > 0) {
      level_external_indices <- points_in_this_level[order(dbscan_result$cluster)]
      level_internal_indices <- as.vector(dbscan_result$cluster)
      num_vertices_in_this_level <- max(level_internal_indices)
    } else {
      level_external_indices <- points_in_this_level
      level_internal_indices <- rep(1, num_points_in_this_level)
      num_vertices_in_this_level <- 1
    } 
  } else if (methods == "pam") {
    level_filter_values <- filter_values[points_in_this_level, , drop = FALSE]
    
    # Check if there are enough points for PAM
    if (nrow(level_filter_values) >= 2) {
      # Choose the number of clusters (k) for PAM, ensure it's valid for the data size
      num_clusters <- min(num_clusters, nrow(level_filter_values) - 1)
      pam_result <- cluster::pam(level_filter_values, k = num_clusters)
      
      # If PAM finds clusters
      if (max(pam_result$clustering) > 0) {
        level_external_indices <- points_in_this_level[order(pam_result$clustering)]
        level_internal_indices <- as.vector(pam_result$clustering)
        num_vertices_in_this_level <- max(level_internal_indices)
      } else {
        level_external_indices <- points_in_this_level
        level_internal_indices <- rep(1, num_points_in_this_level)
        num_vertices_in_this_level <- 1
      }
    } else {
      # If not enough points, treat all as one cluster
      level_external_indices <- points_in_this_level
      level_internal_indices <- rep(1, num_points_in_this_level)
      num_vertices_in_this_level <- 1
    }
  }

  return(list(
    num_vertices = num_vertices_in_this_level,
    external_indices = level_external_indices,
    internal_indices = level_internal_indices
  ))
}

#' Cut the hierarchical clustering tree to define clusters
#'
#' @param heights Heights of the clusters.
#' @param diam Diameter of the clusters.
#' @param num_bins_when_clustering Number of bins when clustering.
#' @return The cutoff height for the clusters.
#' @importFrom graphics hist
#' @export
cluster_cutoff_at_first_empty_bin <- function(heights, diam, num_bins_when_clustering) {
  if (length(heights) == 1) {
    if (heights == diam) {
      cutoff <- Inf
    }
  }
  bin_breaks <- seq(from=min(heights), to=diam, by=(diam - min(heights))/num_bins_when_clustering)
  
  if (length(bin_breaks) == 1) { bin_breaks <- 1 }
  
  myhist <- hist(c(heights,diam), breaks=bin_breaks, plot=FALSE)
  z <- (myhist$counts == 0)
  if (sum(z) == 0) {
    cutoff <- Inf
  } else {
    cutoff <- myhist$mids[ min(which(z == TRUE)) ]
  }
  return(cutoff)
}

#' Find the optimal number of clusters for k-means
#'
#' This function calculates the total within-cluster sum of squares (WSS) for a range
#' of cluster numbers and identifies the best number of clusters (k) based on the 
#' elbow method.
#'
#' @param dist_object A distance matrix or data frame containing the data to be clustered.
#' @param max_clusters The maximum number of clusters to test for k-means. Default is 10.
#' @return The optimal number of clusters (k) based on the elbow method.
#' @importFrom stats kmeans
#' @export
find_best_k_for_kmeans <- function(dist_object, max_clusters = 10) {
  # elbow method
  wss_values <- numeric(max_clusters)
  
  for (k in 1:max_clusters) {
    kmean_result <- kmeans(dist_object, centers = k, nstart = 25)  # nstart for more stable results
    wss_values[k] <- kmean_result$tot.withinss  # Total within-cluster sum of squares
  }
  
  differences <- diff(wss_values)
  second_differences <- diff(differences)
  
  best_k <- which(second_differences == min(second_differences)) + 1
  
  return(best_k)
}