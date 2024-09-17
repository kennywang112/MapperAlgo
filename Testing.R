library(tidyverse)
# library(TDA)
library(networkD3)
library(igraph)
library(ks)
source('R/EdgeVertices.R')
source('R/MapperAlgo.R')

data("iris")

time_taken <- system.time({
  Traffic.mapper <- MapperAlgo(
    filter_values = iris[,1:4],
    intervals = 4,
    percent_overlap = 50,
    num_bins_when_clustering = 30)
})
time_taken

Traffic.graph <- graph.adjacency(Traffic.mapper$adjacency, mode="undirected")
l = length(V(Traffic.graph))
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# Distribution of specific variable in each vertex Majority vote
var.maj.vertex <- c()
filter.vertex <- c()

for (i in 1:l){
  points.in.vertex <- Traffic.mapper$points_in_vertex[[i]]
  Mode.in.vertex <- Mode(iris$Species[points.in.vertex])
  var.maj.vertex <- c(var.maj.vertex,as.character(Mode.in.vertex))
  # filter.vertex <- c(filter.vertex,mean(filter.kde[points.in.vertex]))
}
# Size
vertex.size <- rep(0,l)
for (i in 1:l){
  points.in.vertex <- Traffic.mapper$points_in_vertex[[i]]
  vertex.size[i] <- length((Traffic.mapper$points_in_vertex[[i]]))
}
MapperNodes <- mapperVertices(Traffic.mapper, 1:nrow(iris))
MapperNodes$var.maj.vertex <- as.factor(var.maj.vertex)
# MapperNodes$filter.kde <- filter.vertex
MapperNodes$Nodesize <- vertex.size
MapperLinks <- mapperEdges(Traffic.mapper)
forceNetwork(Nodes = MapperNodes, Links = MapperLinks, Target = "Linktarget",
             Value = "Linkvalue", NodeID = "Nodename", Nodesize = "Nodesize",
             Group = "var.maj.vertex", opacity = 1, zoom = TRUE,
             linkDistance = 10, charge = -10, legend = TRUE)

