% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EvaluatePerformance.R
\name{EvaluatePerformancePerGene}
\alias{EvaluatePerformancePerGene}
\title{Evaluate Performance Per Gene}
\usage{
EvaluatePerformancePerGene(
  outs.metadata,
  spatial,
  multiome,
  spatial_assay = "Xenium",
  multiome_assay = "RNA",
  plot_dir = NULL
)
}
\arguments{
\item{outs.metadata}{A output from `run_projection()`, containing spatial coordinate of multiome cells}

\item{spatial}{A spatial Seurat object}

\item{multiome}{A multiome Seurat object}

\item{spatial_assay}{assay within spatial object shared with multiome data to integrate}

\item{multiome_assay}{assay within multiome object shared with spatial data to integrate}

\item{plot_dir}{Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')}
}
\value{
A dataframe showing Spearman and Pearson correlation between spatial-multiome cells for each gene.
}
\description{
Evaluate the accuracy of integration between spatial and multiome data
by comparing the expression of each gene between spatial-multiome cell counterparts.
}
