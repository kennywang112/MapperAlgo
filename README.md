# sTopological Data Analysis: Mapper Algorithm

This package is based on the `TDAmapper` package by Paul Pearson. You can view the original package [here](https://github.com/paultpearson/TDAmapper). Since the original package hasn't been updated in over seven years, this version is focused on optimization. By incorporating vector computation into the Mapper algorithm, this package aims to significantly improve its performance.

## Get started quickly

![Mapper](man/figures/mapper.png) Step visualize from [Skaf et al.](https://www.sciencedirect.com/science/article/pii/S1532046422000983)

**Mapper is basically a tree-step process:**

1\. **Cover**: This step splits the data into overlapping intervals and creates a cover for the data.

2\. **Cluster**: This step clusters the data points in each interval the cover creates.

3\. **Simplicial Complex**: This step combines the two steps above, which connects the data points in the cover to create a simplicial complex.

> you can know more about the basic here: Chazal, F., & Michel, B. (2021). An introduction to topological data analysis: fundamental and practical aspects for data scientists. Frontiers in artificial intelligence, 4, 667963.

**Besides to the steps above, you can find the following code in the package:**

1.  Mapper.R: Combining the three steps above
2.  ConvertLevelset.R: Converting a Flat Index to a Multi-index, or vice versa.
3.  EdgeVertices.R This is to find the nodes for plot, not for the Mapper algorithm.

## Goals

Although this project serves as a personal training exercise, I have set several key objectives:

1.  **Optimization**: While the current version speeds up computations by 100 times as the dataset grows, there are still some computational challenges that need to be addressed.

2.  **Expanded Clustering Methods**: Clustering is a crucial component of the Mapper algorithm. In addition to hierarchical clustering, I aim to include a variety of clustering techniques to increase flexibility and adaptability.

3.  **Code Structure**: The code is still under development and may be challenging to understand in its current form. My goal is to streamline the structure and provide a simple white paper that explains how to use the method effectively.

## Why R?

While many TDA methods in Python have fewer computational limitations, there is still significant room for improvement in R, especially in terms of optimization and scalability.

## Stay Updated

I've written some articles on Medium, which you can find [here](https://medium.com/@kennywang2003) to get familiar with topological data analysis. I'll be continuously updating my work, and I welcome any feedback!

### Build And Submit:

This is for the author to submit the package to CRAN.

``` r
devtools::build()
devtools::submit_cran()
```
