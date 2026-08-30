#' Perform clustering within a level set using an mlr3cluster Learner
#'
#' @param original_data Original dataframe, not the filter values.
#' @param filter_values The filter values.
#' @param points_in_this_level Points in the current level set.
#' @param method An mlr3cluster Learner, e.g. mlr3cluster::lrn("clust.kmeans", centers = 3).
#' @return A list containing the number of vertices, external indices, and internal indices.
#' @export
perform_clustering_mlr3 <- function(
    original_data, filter_values, points_in_this_level, method
) {
  num_points_in_this_level <- length(points_in_this_level)

  if (num_points_in_this_level == 0) {
    return(list(num_vertices = 0, external_indices = NULL, internal_indices = NULL))
  }

  if (num_points_in_this_level == 1) {
    return(list(num_vertices = 1, external_indices = points_in_this_level, internal_indices = c(1)))
  }

  if (inherits(method, "LearnerClustAgnes")) {
    stop("Hierarchical clustering is not yet supported via `method`; use `methods = \"hierarchical\"` instead.")
  }

  level_data <- original_data[points_in_this_level, , drop = FALSE]

  cluster_vector <- tryCatch({
    task <- mlr3cluster::TaskClust$new(id = "level", backend = level_data)
    learner <- method$clone(deep = TRUE)
    learner$train(task)
    as.vector(learner$predict(task)$partition)
  }, error = function(e) {
    rep(1, num_points_in_this_level)
  })

  list(
    num_vertices = max(cluster_vector),
    external_indices = points_in_this_level,
    internal_indices = cluster_vector
  )
}
