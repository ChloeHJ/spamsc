% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConvertintoSpatialObject.R
\name{ConvertAnndata}
\alias{ConvertAnndata}
\title{Convert Seurat to Anndata}
\usage{
ConvertAnndata(
  seurat,
  assays = c("RNA"),
  dimreducs = c("umap"),
  filename = "spamscOutput"
)
}
\arguments{
\item{seurat}{A Seurat object. Possibly output from `MultiometoSpatialObject()`}

\item{assays}{An assay to convert to Anndata}

\item{dimreducs}{Only keep a subset of DimReducs specified here (if NULL, remove all DimReducs)}

\item{filename}{name of the file to save Anndata}
}
\value{
A h5ad Anndata
}
\description{
Processing function to convert a Seurat object into an Anndata
}
