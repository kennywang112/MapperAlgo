#' Construct adjacency matrix of the simplicial complex
#'
#' @param filter_values A matrix of filter values.
#' @param vertex_index The number of vertices.
#' @param num_levelsets The total number of level sets.
#' @param num_intervals A vector representing the number of intervals for each filter.
#' @param vertices_in_level_set A list where each element contains the vertices corresponding to each level set.
#' @param points_in_vertex A list where each element contains the points corresponding to each vertex.
#' @return An adjacency matrix representing the simplicial complex.
#' @export
simplcial_complex <- function(
    filter_values, vertex_index, num_levelsets, num_intervals, vertices_in_level_set, points_in_vertex
) {
  filter_output_dim <- dim(filter_values)[2] # columns
  # create empty adjacency matrix to store the connections between vertices
  adja <- mat.or.vec(vertex_index, vertex_index)
  for (lsfi in 1:num_levelsets) {

    lsmi <- to_lsmi(lsfi, num_intervals)
    # Find adjacent level sets +1 of each entry in lsmi (within bounds of num_intervals)
    # Need to_lsfi to do this easily.
    for (k in 1:filter_output_dim) {
      # check admissibility condition
      if (lsmi[k] >= num_intervals[k]) { next }
      lsmi_adjacent <- lsmi + diag(filter_output_dim)[, k]
      lsfi_adjacent <- to_lsfi(lsmi_adjacent, num_intervals)
      
      v1_set <- vertices_in_level_set[[lsfi]]
      v2_set <- vertices_in_level_set[[lsfi_adjacent]]
      
      if (length(v1_set) < 1 | length(v2_set) < 1) { next }
      # construct adjacency matrix
      for (v1 in v1_set) {
        for (v2 in v2_set) {
          adja[v1, v2] <- (length(intersect(
            points_in_vertex[[v1]], points_in_vertex[[v2]])) > 0)
          
          adja[v2, v1] <- adja[v1,v2]
        }
      }
    }
  }
  return(adja)
}