% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ComputeSpatialNiche.R
\name{FindNiches}
\alias{FindNiches}
\title{FindNiches}
\usage{
FindNiches(
  spatial,
  niche.matrix,
  resolution = 0.1,
  dims = 1:10,
  niche.assay.name = "niche",
  harmonise = F,
  harmony.group.by = NULL,
  niche.name = "TissueNiche"
)
}
\arguments{
\item{spatial}{seurat spatial}

\item{niche.matrix}{any niche matrix that matrches the cells of the seurat}

\item{resolution}{- clustering resolution, this generally should be lower than for single cell}

\item{dims}{- how many dims to use, again normally 10 is good, too many gets weird small clusters}

\item{niche.assay.name}{name of the new assay where teh data will be stored in the seurat}

\item{harmonise}{whether to harmonise the niche matrix between any variables -e.g. different slides}

\item{harmony.group.by}{if harmonise, provide group.by arguement}

\item{niche.name}{what name to give to the niche groups in the spatial meta data}
}
\description{
Find niche clusters given an arbitraty niche matrix - can be expression or cell type counts matrix
}
