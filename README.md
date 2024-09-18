# Topological Data Analysis: Mapper Algorithm

This package is based on the `TDAmapper` package by Paul Pearson. You can view the original package [here](https://github.com/paultpearson/TDAmapper). Since the original package hasn't been updated in over seven years, this version is focused on optimization. By incorporating vector computation into the Mapper algorithm, this package aims to significantly improve its performance.

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
