#' Create Mapper Edges
#'
#' This function generates the edges of the Mapper graph by analyzing the adjacency matrix.
#' It returns a data frame with source and target vertices that are connected by edges.
#'
#' @param m The Mapper output object that contains the adjacency matrix and other graph components.
#' @return A data frame containing the source (`Linksource`), target (`Linktarget`), and edge values (`Linkvalue`) for the graph's edges.
#' @export
mapperEdges <- function(m) {
  linksource <- c()
  linktarget <- c()
  linkvalue <- c()
  k <- 1
  for (i in 2:m$num_vertices) {
    for (j in 1:(i-1)) {
      if (m$adjacency[i,j] == 1) {
        linksource[k] <- i - 1
        linktarget[k] <- j - 1
        linkvalue[k] <- 2
        k <- k + 1
      }
    }
  }
  return(data.frame(Linksource = linksource,
                    Linktarget = linktarget, 
                    Linkvalue = linkvalue))
}

#' Create Mapper Vertices
#'
#' This function generates the vertices of the Mapper graph, including their labels and groupings.
#' It returns a data frame with the vertex names, the group each vertex belongs to, and the size of each vertex.
#'
#' @param m The Mapper output object that contains information about the vertices and level sets.
#' @param pt_labels A vector of point labels to be assigned to the points in each vertex.
#' @return A data frame containing the vertex names (`Nodename`), group information (`Nodegroup`), and vertex sizes (`Nodesize`).
#' @export
mapperVertices <- function(m, pt_labels) {
  labels_in_vertex <- lapply(m$points_in_vertex, FUN = function(v) { pt_labels[v] })
  nodename <- sapply(sapply(labels_in_vertex, as.character), paste0, collapse = ", ")
  nodename <- paste0("V", 1:m$num_vertices, ": ", nodename)
  
  nodegroup <- m$level_of_vertex
  nodesize <- sapply(m$points_in_vertex, length)
  
  return(data.frame(Nodename = nodename, 
                    Nodegroup = nodegroup, 
                    Nodesize = nodesize))
}
