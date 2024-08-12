#' @title Compute coexpression
#' @description Compute correlation across genes of interest
#'
#' @param data gene x cell matrix containing expression values
#'
#' @return A gene x gene matrix showing Pearson correlation coefficient
#' @export
ComputeCoexpression <- function(data ){

  cor_mtx <- cor(data)
  #heatmap(cor_mtx, scale = "row")

  return(cor_mtx)
}

#' @title Visualise coexpression
#' @description Plot heatmap showing the correlation between genes.
#'
#' @param cor_mtx output correlation matrix from `compute_coexpression()`
#' @param cor.threshold A minimum Pearson correlation coefficient threshold to plot genes. Default: 0.4
#'
#' @return A heatmap showing Perason correalaiton coefficient across genes
#' @export
VisualiseCoexpression <- function(cor_mtx,  cor.threshold = 0.4){

  # Set the correlation coefficient threshold
  high_cor_genes <- which(abs(cor_mtx) > cor.threshold & upper.tri(cor_mtx, diag = FALSE), arr.ind = TRUE)
  high_cor_gene_names <- rownames(cor_mtx)[high_cor_genes[, 1]]
  high_cor_gene_names <- unique(c(high_cor_gene_names, rownames(cor_mtx)[high_cor_genes[, 2]]))
  print('Genes co-expressed beyond threshold:')
  print(high_cor_gene_names)

  print(pheatmap(cor_mtx[high_cor_gene_names, high_cor_gene_names], cluster_rows = TRUE, cluster_cols = TRUE))

  plot_mtx <- cor_mtx[high_cor_gene_names, high_cor_gene_names]
  return(plot_mtx)

  return()
}


#' @title Visualise coexpression within range
#' @description Plot heatmap showing the correlation between genes.
#'
#' @param cor_mtx output correlation matrix from `compute_coexpression()`
#' @param cor.top.threshold A maximum Pearson correlation coefficient threshold to plot genes. Default: 0.4
#' @param cor.bottom.threshold A minimum Pearson correlation coefficient threshold to plot genes. Default: 0.4
#' @return A heatmap showing Perason correalaiton coefficient across genes
#'
#' @export
VisualiseCoexpressionWindow <- function(cor_mtx,  cor.top.threshold = 0.8, cor.bottom.threshold = 0.4){

  # Find high correlation genes
  high_cor_genes <- which(abs(cor_mtx) > cor.bottom.threshold & abs(cor_mtx) < cor.top.threshold & upper.tri(cor_mtx, diag = FALSE), arr.ind = TRUE)

  # Extract gene names
  high_cor_gene_names <- unique(c(rownames(cor_mtx)[high_cor_genes[, 1]], rownames(cor_mtx)[high_cor_genes[, 2]]))

  # Print genes co-expressed within threshold
  print('Genes co-expressed within threshold:')
  print(high_cor_gene_names)

  # Plot heatmap for filtered genes
  if (length(high_cor_gene_names) > 1) {
    plot_mtx <- cor_mtx[high_cor_gene_names, high_cor_gene_names]
    pheatmap(plot_mtx, cluster_rows = TRUE, cluster_cols = TRUE)
  } else {
    print('Not enough genes to plot a heatmap.')
  }

  return(plot_mtx)


}

