library(networkD3)
library(igraph)

source('R/EdgeVertices.R')
source('R/ConvertLevelsets.R')
source('R/Cover.R')
source('R/Cluster.R')
source('R/SimplicialComplex.R')
source('R/MapperAlgo.R')

data("iris")

time_taken <- system.time({
  Mapper <- MapperAlgo(
    filter_values = iris[,1:4],
    intervals = 4,
    percent_overlap = 50,
    num_bins_when_clustering = 30,
    methods ="kmeans")
})
time_taken

Graph <- graph.adjacency(Mapper$adjacency, mode="undirected")
l = length(V(Graph))
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# Distribution of specific variable in each vertex Majority vote
var.maj.vertex <- c()
filter.vertex <- c()

for (i in 1:l){
  points.in.vertex <- Mapper$points_in_vertex[[i]]
  Mode.in.vertex <- Mode(iris$Species[points.in.vertex])
  var.maj.vertex <- c(var.maj.vertex,as.character(Mode.in.vertex))
  # filter.vertex <- c(filter.vertex,mean(filter.kde[points.in.vertex]))
}
# Size
vertex.size <- rep(0,l)
for (i in 1:l){
  points.in.vertex <- Mapper$points_in_vertex[[i]]
  vertex.size[i] <- length((Mapper$points_in_vertex[[i]]))
}
MapperNodes <- mapperVertices(Mapper, 1:nrow(iris))
MapperNodes$var.maj.vertex <- as.factor(var.maj.vertex)
# MapperNodes$filter.kde <- filter.vertex
MapperNodes$Nodesize <- vertex.size
MapperLinks <- mapperEdges(Mapper)
forceNetwork(Nodes = MapperNodes, Links = MapperLinks, Target = "Linktarget",
             Value = "Linkvalue", NodeID = "Nodename", Nodesize = "Nodesize",
             Group = "var.maj.vertex", opacity = 1, zoom = TRUE,
             linkDistance = 10, charge = -10, legend = TRUE)

