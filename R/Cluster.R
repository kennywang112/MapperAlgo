#' Perform clustering within a level set
#'
#' @param points_in_this_level Points in the current level set.
#' @param filter_values The filter values.
#' @param num_bins_when_clustering Number of bins when clustering.
#' @return A list containing the number of vertices, external indices, and internal indices.
#' @importFrom stats as.dist hclust cutree dist
#' @export
perform_clustering <- function(points_in_this_level, filter_values, num_bins_when_clustering) {
  num_points_in_this_level <- length(points_in_this_level)
  
  if (num_points_in_this_level == 0) {
    return(list(num_vertices = 0, external_indices = NULL, internal_indices = NULL))
  }
  
  if (num_points_in_this_level == 1) {
    return(list(num_vertices = 1, external_indices = points_in_this_level, internal_indices = c(1)))
  }
  
  level_dist_object <- as.dist(as.matrix(dist(filter_values))[points_in_this_level, points_in_this_level])
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
