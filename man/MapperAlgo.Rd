\name{MapperAlgo}
\alias{MapperAlgo}
\title{
Topological data analysis: Mapper algorithm
}
\description{
The Mapper algorithm is a method for topological data analysis that provides a way to visualize the structure of high-dimensional data. The Mapper algorithm is a generalization of the Reeb graph construction, which is a method for visualizing the topology of scalar fields.
}
\usage{
MapperAlgo(filter_values, intervals, percent_overlap, num_bins_when_clustering)
}
\arguments{
  \item{filter_values}{A data frame or matrix of the data to be analyzed.}
  \item{intervals}{An integer specifying the number of intervals to divide the filter values into.}
  \item{percent_overlap}{An integer specifying the percentage of overlap between consecutive intervals.}
  \item{num_bins_when_clustering}{An integer specifying the number of bins to use when clustering the data.}
}
\value{
  An adjacency matrix and other components of the Mapper graph, including:
  \item{adjacency}{An adjacency matrix of the Mapper graph.}
  \item{num_vertices}{The number of vertices in the Mapper graph.}
  \item{level_of_vertex}{A vector specifying the level of each vertex.}
  \item{points_in_vertex}{A list of the indices of the points in each vertex.}
  \item{points_in_level_set}{A list of the indices of the points in each level set.}
  \item{vertices_in_level_set}{A list of the indices of the vertices in each level set.}
}
\references{
The original paper on the Mapper algorithm is:
G. Singh, F. Memoli, G. Carlsson (2007). Topological Methods for the Analysis of High Dimensional Data Sets and 3D Object Recognition, Point Based Graphics 2007, Prague, September 2007.
This code is based on Paul Pearson's implementation of the Mapper algorithm in R, optimized for speed and memory usage. You can install using the following command: devtools::install_github("paultpearson/TDAmapper")
}
\author{
ChiChien Wang
}
\examples{
library(igraph)

data("iris")

mapper <- MapperAlgo(
  filter_values = iris[,1:4],
  intervals = 4,
  percent_overlap = 50,
  num_bins_when_clustering = 30)
    
graph <- graph.adjacency(mapper$adjacency, mode="undirected")
l = length(V(graph))
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# Distribution of specific variable in each vertex - Majority vote
var.maj.vertex <- c()
filter.vertex <- c()

for (i in 1:l){
  points.in.vertex <- mapper$points_in_vertex[[i]]
  Mode.in.vertex <- Mode(iris$Species[points.in.vertex])
  var.maj.vertex <- c(var.maj.vertex, as.character(Mode.in.vertex))
}

# Size
vertex.size <- rep(0, l)
for (i in 1:l){
  points.in.vertex <- mapper$points_in_vertex[[i]]
  vertex.size[i] <- length(mapper$points_in_vertex[[i]])
}

MapperNodes <- mapperVertices(mapper, 1:nrow(iris))
MapperNodes$var.maj.vertex <- as.factor(var.maj.vertex)
MapperNodes$Nodesize <- vertex.size
MapperLinks <- mapperEdges(mapper)
}
\keyword{mapper}
