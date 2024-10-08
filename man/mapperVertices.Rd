% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EdgeVertices.R
\name{mapperVertices}
\alias{mapperVertices}
\title{Create Mapper Vertices}
\usage{
mapperVertices(m, pt_labels)
}
\arguments{
\item{m}{The Mapper output object that contains information about the vertices and level sets.}

\item{pt_labels}{A vector of point labels to be assigned to the points in each vertex.}
}
\value{
A data frame containing the vertex names (\code{Nodename}), group information (\code{Nodegroup}), and vertex sizes (\code{Nodesize}).
}
\description{
This function generates the vertices of the Mapper graph, including their labels and groupings.
It returns a data frame with the vertex names, the group each vertex belongs to, and the size of each vertex.
}
