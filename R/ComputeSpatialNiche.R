#' @title BuildNicheExpressionAssay
#' @description Creates a matrix of expression values for each cell's nearest spatial neighbours
#' To be used as an alternative niche assay without needing to specify cell types apriori - can test for enrichment afterwards
#' So this is a more unsupervised approach
#' @param coords cell x 2 matrix containing spatial coordinate of cells
#' @param data gene x cell matrix containing expression values
#'@param neighbors.k nearest neighbours to use
#'@return matrix of nearest neighbour aggregated counts for each feature, excluding the cell itself
#' @export
BuildNicheExpressionAssay <- function(
    data,
    coords,
    neighbors.k = 30
) {

  cells <- colnames(data)
  neighbors <- FindNeighbors(coords, k.param = neighbors.k)
  neighbors$nn <- neighbors$nn[cells, cells]

  diag(neighbors$nn) <- 0 # dont count transcriptome of the cell itself, just neighbours?

  mt <- data[, cells]

  sum.mtx <- as.matrix(neighbors$nn %*% t(mt))

  return(sum.mtx)
}



#' @title FindNiches
#' @description Find niche clusters given an arbitraty niche matrix - can be expression or cell type counts matrix
#' @param spatial seurat spatial
#' @param niche.matrix any niche matrix that matrches the cells of the seurat
#' @param resolution - clustering resolution, this generally should be lower than for single cell
#' @param dims - how many dims to use, again normally 10 is good, too many gets weird small clusters
#' @param niche.assay.name name of the new assay where teh data will be stored in the seurat
#' @param harmonise whether to harmonise the niche matrix between any variables -e.g. different slides
#' @param harmony.group.by if harmonise, provide group.by arguement
#' @param niche.name what name to give to the niche groups in the spatial meta data
#' @export
FindNiches <- function(spatial, niche.matrix, resolution=0.1,
                       dims=1:10, niche.assay.name="niche",
                       harmonise=F, harmony.group.by=NULL, niche.name="TissueNiche"){


  spatial[[niche.assay.name]] <- CreateAssayObject(niche.matrix[, Cells(spatial)])
  spatial@active.assay <- niche.assay.name
  spatial <- SCTransform(spatial, assay = niche.assay.name, clip.range = c(-10, 10))
  # seq.gs <- NormalizeData(spatial, assay = niche.assay.name)
  # seq.gs <- ScaleData(spatial, assay = niche.assay.name)
  # seq.gs <- FindVariableFeatures(spatial, assay = niche.assay.name)
  spatial <- RunPCA(spatial, npcs = 30)
  tmp <- spatial@active.ident

  if( harmonise){

    spatial <- harmony::RunHarmony(spatial, group.by.vars=harmony.group.by, assay=niche.assay.name)
    spatial <- FindNeighbors(spatial, reduction = "harmony", dims = dims)
    spatial <- FindClusters(spatial, resolution = resolution)
  }
  else{

    spatial <- FindNeighbors(spatial, reduction = "pca", dims = dims)
    spatial <- FindClusters(spatial, resolution = resolution)
  }

  spatial <- StashIdent(spatial, niche.name) # stack the niche cluster identity in the meta data
  spatial@active.ident <- tmp ## reset back to original active ident so it is not accidentally lost, like I am likely to do
  return(spatial)
}

#' @title CellTypeNicheEnrichment
#' @description Given a vector of niches and a corresponding vector of cell types, will calculate cell type enrichment table in each niche
#' to be used with downstream visualisation functions
#' @param niches - vector of niche names per cell
#' @param cell.types - vector of cell types per cell, in the same cell order as that of niches
#' @return data frame with enrichment estimate and p value for each  cell type - cell niche combination
#' @export
CellTypeNicheEnrichment <- function(niches, cell.types){

  df <- list()

  for ( niche in unique(niches)){

    for( cell in unique(cell.types)){

      res <- fisher.test(x=niches == niche, y=cell.types == cell )

      df[[paste0(cell, niche)]] <- c(res$estimate, res$p.value, cell, niche)
    }
  }

  df <- as.data.frame(do.call(rbind, df))
  colnames(df) <- c("Estimate", "p.val", "Cell Type", "Niche")
  df$p.val <- as.numeric(df$p.val)
  df$FDR <- p.adjust(df$p.val)

  return(df)
}



#' @title PlotNicheEnrichmentBarPlot
#' @description given niche results from CellTypeNicheEnrichment, will plot a cell type enrichment barplot for any query niche
#' @param niche.outs - output table from `CellTypeNicheEnrichment`
#' @param niche - which niche to plot, must be present in the results
#' @return barplot with enrichment results
#'
#' @export
PlotNicheEnrichmentBarPlot <- function(niche.outs, niche=0){

  df <- niche.outs[niche.outs$Niche == niche,]
  df <- df[order(as.numeric(df$Estimate)), ]
  df$`Cell Type` <- factor(df$`Cell Type`, levels=df$`Cell Type`)
  ggplot(df, aes(`Cell Type`, log2(as.numeric(Estimate)),
                 fill=-log10(as.numeric(FDR) + 10^-300))
  ) + geom_bar(stat="identity", colour="black"
  )  + coord_flip() + scale_fill_viridis_c() + theme_classic(base_size = 16) + labs(
    x="Cell Type", y="Relative Depletion <-----> Relative Enrichment", fill="-log10 FDR"
  ) + geom_hline(yintercept = c(-1, 1), lty=2, colour="red")
}

#' @title PlotCellTypeEnrichmentBarPlot
#' @description Given niche results from CellTypeNicheEnrichment, will plot a niche enrichment barplot for any query cell type
#' @param niche.outs - output table from `CellTypeNicheEnrichment`
#' @param CellType - which cell type to plot, must be present in the results
#' @return barplot with enrichment results
#'
#' @export
PlotCellTypeEnrichmentBarPlot <- function(niche.outs, CellType=0){

  df <- niche.outs[niche.outs$`Cell Type` == CellType,]
  df <- df[order(as.numeric(df$Estimate)), ]
  df$Niche <- factor(df$Niche, levels=sort(as.numeric(as.character(unique(df$Niche)))))
  ggplot(df, aes(Niche, log2(as.numeric(Estimate)),
                 fill=-log10(as.numeric(FDR) + 10^-300))
  ) + geom_bar(stat="identity", colour="black"
  )  + coord_flip() + scale_fill_viridis_c() + theme_classic(base_size = 16) + labs(
    x="Niche", y="Relative Depletion <-----> Relative Enrichment", fill="-log10 FDR"
  ) + geom_hline(yintercept = c(-1, 1), lty=2, colour="red")
}



#' @title PlotNicheEnrichmentHeatmap
#' @description Given niche results from cellTypeNicheEnrichment, will plot an overview heatmap of cell type to cell niche enrichment.
#' @param niche.outs  - output table from `CellTypeNicheEnrichment` function
#' @return heatmap with enrichment results
#'
#' @export
PlotNicheEnrichmentHeatmap <- function(niche.outs){

  niche.outs$Estimate <- as.numeric(niche.outs$Estimate)
  mat <- log2(reshape2::acast(niche.outs, `Cell Type`~ Niche, value.var = "Estimate")+ 0.01)
  niche.outs$STAR <- "n.s"; niche.outs$STAR[niche.outs$FDR < 0.05] <- "*"; niche.outs$STAR[niche.outs$FDR < 0.01] <- "**"; niche.outs$STAR[niche.outs$FDR < 0.001] <- "***"
  mat.labels <-reshape2::acast(niche.outs, `Cell Type`~ Niche, value.var = "STAR")
  mat[mat < -4] <- -4; mat[ mat > 4] <- 4
  pheatmap::pheatmap(mat, border_color = "black", display_numbers = mat.labels, fontsize_number = 14,
                     color=SpatialColors(n=100), number_color = "black" )

}

