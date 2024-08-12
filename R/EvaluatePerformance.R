#' @title Evaluate Performance Per Gene
#' @description Evaluate the accuracy of integration between spatial and multiome data
#' by comparing the expression of each gene between spatial-multiome cell counterparts.
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param spatial_assay assay within spatial object shared with multiome data to integrate
#' @param multiome_assay assay within multiome object shared with spatial data to integrate
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return A dataframe showing Spearman and Pearson correlation between spatial-multiome cells for each gene.
#' @export
EvaluatePerformancePerGene <- function(outs.metadata, spatial, multiome, spatial_assay = 'Xenium', multiome_assay = 'RNA', plot_dir = NULL){

  DefaultAssay(spatial) <- spatial_assay
  DefaultAssay(multiome) <- multiome_assay

  id_match <- outs.metadata %>% dplyr::select(multiome_id, spatial_id)
  print(paste0('# of mulitome cells without spatial counterpart: ', length(which(is.na(id_match$spatial_id))) ))
  ids <- id_match %>% filter(!is.na(spatial_id)) # those with spatial counterparts

  # scale data on all features
  multiome <- ScaleData(multiome, features = rownames(multiome))
  spatial <- ScaleData(spatial, features = rownames(spatial))

  # for each gene - each point represents a gene
  cor_pr_gene <- list(); cor_sp_gene <- list()
  shared_genes <- intersect(rownames(spatial@assays[[spatial_assay]]@scale.data), rownames(multiome@assays[[multiome_assay]]@scale.data))
  for(g in shared_genes){
    multiome_exp <- as.matrix(multiome@assays[[multiome_assay]]@scale.data[g, ids$multiome_id ])
    spatial_exp <- as.matrix(spatial@assays[[spatial_assay]]@scale.data[g, ids$spatial_id ])
    cor_pr_gene[g] <- cor(multiome_exp, spatial_exp, method = 'pearson', use = "complete.obs")
    cor_sp_gene[g] <- cor(multiome_exp, spatial_exp, method = 'spearman', use = "complete.obs")
  }

  cor_gene <- cbind(t(data.frame(cor_pr_gene)), t(data.frame(cor_sp_gene)))
  colnames(cor_gene) <- c('pearson', 'spearman')

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.histogram.per.gene.pearson.cor.pdf'), width = 4, height = 3.5)}
  print(ggplot(data = as.data.frame(cor_gene), aes(x = pearson)) + geom_histogram() + theme_classic() + labs(x = 'Pearson correlation', y = '# of Genes') +
          ggtitle('Histogram of per gene\npearson correlation coefficient') + theme(plot.title = element_text(hjust = 0.5)))
  if(!is.null(plot_dir)){dev.off()}

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.histogram.per.gene.spearman.cor.pdf'), width = 4, height = 3.5)}
  print(ggplot(data = as.data.frame(cor_gene), aes(x = spearman)) + geom_histogram() + theme_classic() + labs(x = 'Spearman correlation', y = '# of Genes') +
          ggtitle('Histogram of per gene\nspearman correlation coefficient') + theme(plot.title = element_text(hjust = 0.5)))
  if(!is.null(plot_dir)){dev.off()}

  print('Summary of per gene spearman correlation coefficient')
  print(summary(cor_gene[,2]) )


  outs <- cor_gene
  return(outs)

}


#' @title Evaluate Performance Per Cell
#' @description Evaluate the accuracy of integration between spatial and multiome data
#' by comparing the expression of spatial-multiome cell counterparts,
#' and compare to randomly shuffled background
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param spatial_assay assay within spatial object shared with multiome data to integrate
#' @param multiome_assay assay within multiome object shared with spatial data to integrate
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return A dataframe showing Spearman and Pearson correlation between spatial-multiome cells
#' @export
EvaluatePerformancePerCell <- function(outs.metadata, spatial, multiome, spatial_assay = 'Xenium', multiome_assay = 'RNA', plot_dir = NULL){

  DefaultAssay(spatial) <- spatial_assay
  DefaultAssay(multiome) <- multiome_assay

  id_match <- outs.metadata %>% dplyr::select(multiome_id, spatial_id)
  print(paste0('# of mulitome cells without spatial counterpart: ', length(which(is.na(id_match$spatial_id))) ))
  ids <- id_match %>% filter(!is.na(spatial_id)) # those with spatial counterparts

  # scale data on all features
  multiome <- ScaleData(multiome, features = rownames(multiome))
  spatial <- ScaleData(spatial, features = rownames(spatial))
  shared_genes <- intersect(rownames(spatial), rownames(multiome))

  # for each cell - each point represents a cell
  cor_pr_cell <- list(); cor_sp_cell <- list()
  for(i in 1:nrow(ids)){
    multiome_exp <- as.matrix(multiome@assays[[multiome_assay]]@scale.data[shared_genes, ids$multiome_id[i] ])
    spatial_exp <- as.matrix(spatial@assays[[spatial_assay]]@scale.data[shared_genes,ids$spatial_id[i] ])
    cor_pr_cell[i] <- cor(multiome_exp, spatial_exp, method = 'pearson', use = "complete.obs")
    cor_sp_cell[i] <- cor(multiome_exp, spatial_exp, method = 'spearman', use = "complete.obs")
  }

  cor_cell <- cbind(t(data.frame(cor_pr_cell)), t(data.frame(cor_sp_cell)))
  colnames(cor_cell) <- c('pearson', 'spearman')

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.histogram.per.cell.pearson.cor.pdf'), width = 4, height = 3.5)}
  print(ggplot(data = as.data.frame(cor_cell), aes(x = pearson)) + geom_histogram() + theme_classic() + labs(x = 'Pearson correlation', y = '# of Cells') +
          ggtitle('Histogram of per cell\npearson correlation coefficient') + theme(plot.title = element_text(hjust = 0.5)))
  if(!is.null(plot_dir)){dev.off()}

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.histogram.per.cell.spearman.cor.pdf'), width = 4, height = 3.5)}
  print(ggplot(data = as.data.frame(cor_cell), aes(x = spearman)) + geom_histogram() + theme_classic() + labs(x = 'Spearman correlation', y = '# of Cells') +
          ggtitle('Histogram of per cell\nspearman correlation coefficient') + theme(plot.title = element_text(hjust = 0.5)))
  if(!is.null(plot_dir)){dev.off()}

  print('Summary of per cell spearman correlation coefficient')
  print(summary(cor_cell[,2]) )

  # compare to background (whole spatial slide)
  spatial <- ScaleData(spatial, features = rownames(spatial))
  set.seed(2345);  if(length(colnames(spatial)) < nrow(ids)){randomised_spatial_id <- sample(colnames(spatial), size = nrow(ids), replace = TRUE)
  }else{randomised_spatial_id <- sample(colnames(spatial), size = nrow(ids))}
  cor_sp_random_cell <- list()
  for(i in 1:nrow(ids)){
    multiome_exp <- as.matrix(multiome@assays[[multiome_assay]]@scale.data[shared_genes, ids$multiome_id[i] ])
    spatial_exp <- as.matrix(spatial@assays[[spatial_assay]]@scale.data[shared_genes,randomised_spatial_id[i] ])
    cor_sp_random_cell[i] <- cor(multiome_exp, spatial_exp, method = 'spearman', use = "complete.obs")
  }

  cor_random_cell <- cbind(t(data.frame(cor_sp_cell)), t(data.frame(cor_sp_random_cell)))
  colnames(cor_random_cell) <- c('spearman', 'spearman_bg')

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.histogram.per.cell.spearman.cor.bg.pdf'), width = 5, height = 3.5)}
  print(ggplot(data = as.data.frame(cor_random_cell) %>% melt(), aes(x = value, color = variable)) +
          geom_histogram(fill="white", alpha=0.5, position="identity") +
          theme_classic() + labs(x = 'Spearman correlation', y = '# of Cells') +
          ggtitle('Histogram of per cell spearman correlation coefficient,\ncompared to background') + theme(plot.title = element_text(hjust = 0.5)))
  if(!is.null(plot_dir)){dev.off()}

  print('Summary of background spearman correlation coefficient')
  print(summary(cor_random_cell[, 2]))


  outs <- list( per_cell_correlation = cor_cell,
               per_cell_correlation_bg = cor_random_cell)

  return(outs)

}

#' @title Plot Gene Correlation
#' @description Plot correlation in gene expression between spatial-multiome cell counterparts.
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param spatial_assay assay within spatial object shared with multiome data to integrate
#' @param multiome_assay assay within multiome object shared with spatial data to integrate
#' @param gene gene of interest
#' @param pt.size size of points in scatter plot
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return A scatterplot showing gene expression in spatial and multiome data
#' @export
PlotGeneCorrelation <- function(outs.metadata, gene = 'CD3D', multiome, spatial, spatial_assay = 'Xenium',
                                multiome_assay = 'RNA', pt.size = 0.3, plot_dir = NULL ){
  id_match <- outs.metadata %>% dplyr::select(multiome_id, spatial_id)
  length(which(is.na(id_match$spatial_id))) #multiome cells don't have spatial counterparts
  ids <- id_match %>% filter(!is.na(spatial_id)) # those with spatial counterparts

  DefaultAssay(spatial) <- spatial_assay
  DefaultAssay(multiome) <- multiome_assay

  # scale data on all features
  multiome <- ScaleData(multiome, features = rownames(multiome))
  spatial <- ScaleData(spatial, features = rownames(spatial))

  # per gene expression
  g <- gene
  multiome_g_exp <- as.matrix(multiome@assays[[multiome_assay]]@scale.data[g, ids$multiome_id ])
  spatial_g_exp <- as.matrix(spatial@assays[[spatial_assay]]@scale.data[g, ids$spatial_id ])

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.scatterplot.per.gene.expression.', gene, '.pdf'), width = 2.5, height = 2.5)}
  print(cbind(multiome_g_exp, spatial_g_exp) %>% as.data.frame() %>%
    ggplot( aes(V1, V2)) + geom_point(size = pt.size)  + theme_classic() +
    labs(x = 'Multiome expression', y = 'Spatial expression', title = g) + stat_cor(method="pearson") )
  if(!is.null(plot_dir)){dev.off()}

}

#' @title Plot correlation on spatial coordinate
#' @description Plot correlation in gene expression between spatial-multiome cell counterparts
#' on their spatial coordinates
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param spatial_assay assay within spatial object shared with multiome data to integrate
#' @param multiome_assay assay within multiome object shared with spatial data to integrate
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return A dataframe showing Spearman and Pearson correlation of each multiome cells
#' @export
PlotCorrelationSpatial <- function(outs.metadata, multiome, spatial, spatial_assay = 'Xenium', multiome_assay = 'RNA', plot_dir = NULL){

  DefaultAssay(spatial) <- spatial_assay
  DefaultAssay(multiome) <- multiome_assay
  shared_genes <- intersect(rownames(spatial), rownames(multiome))

  id_match <- outs.metadata %>% dplyr::select(multiome_id, spatial_id)
  length(which(is.na(id_match$spatial_id))) #multiome cells don't have spatial counterparts
  ids <- id_match %>% filter(!is.na(spatial_id)) # those with spatial counterparts

  # scale data on all features
  multiome <- ScaleData(multiome, features = rownames(multiome))
  spatial <- ScaleData(spatial, features = rownames(spatial))

  # for each cell - each point represents a cell
  cor_pr_cell <- list(); cor_sp_cell <- list()
  for(i in 1:nrow(ids)){
    multiome_exp <- as.matrix(multiome@assays[[multiome_assay]]@scale.data[shared_genes, ids$multiome_id[i] ])
    spatial_exp <- as.matrix(spatial@assays[[spatial_assay]]@scale.data[shared_genes,ids$spatial_id[i] ])
    cor_pr_cell[i] <- cor(multiome_exp, spatial_exp, method = 'pearson', use = "complete.obs")
    cor_sp_cell[i] <- cor(multiome_exp, spatial_exp, method = 'spearman', use = "complete.obs")
  }

  cor_cell <- cbind(t(data.frame(cor_pr_cell)), t(data.frame(cor_sp_cell)))
  colnames(cor_cell) <- c('pearson', 'spearman')

  # show correlation on spatial slide
  cor_cell_ids <- cbind(ids, cor_cell)
  outs <- outs.metadata %>% dplyr::select(multiome_id, x, y) %>% left_join(cor_cell_ids, by = 'multiome_id')

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.spatial.coord.per.cell.pearson.correlation.pdf'), width = 5, height = 4)}
  print(outs %>%
          ggplot( aes(y, x, color = pearson)) + geom_point(size = 0.8)  + theme_classic() + scale_color_viridis_c(option = "magma", direction = -1) +
          ggtitle('Per cell pearson correlation coefficient\non spatial coordinates') + theme(plot.title = element_text(hjust = 0.5)))
  if(!is.null(plot_dir)){dev.off()}

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.spatial.coord.per.cell.spearman.correlation.pdf'), width = 5, height = 4)}
  print(outs %>%
          ggplot( aes(y, x, color = spearman)) + geom_point(size = 0.8) + theme_classic() + scale_color_viridis_c(option = "magma", direction = -1) +
          ggtitle('Per cell spearman correlation coefficient\non spatial coordinates') + theme(plot.title = element_text(hjust = 0.5)))
  if(!is.null(plot_dir)){dev.off()}

  return(outs)

}

#' @title Plot summed expression on spatial coordinate
#' @description Plot summed gene counts for each multiome cells on their spatial coordinates
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param spatial_assay assay within spatial object shared with multiome data to integrate
#' @param multiome_assay assay within multiome object shared with spatial data to integrate
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return A dataframe showing summed gene counts of each multiome cells
#' @export
PlotSummedExpressionSpatial <- function(outs.metadata, multiome, spatial, spatial_assay = 'Xenium', multiome_assay = 'RNA', plot_dir = NULL){

  DefaultAssay(spatial) <- spatial_assay
  DefaultAssay(multiome) <- multiome_assay
  shared_genes <- intersect(rownames(spatial), rownames(multiome))

  id_match <- outs.metadata %>% select(multiome_id, spatial_id)
  length(which(is.na(id_match$spatial_id))) #multiome cells don't have spatial counterparts
  ids <- id_match %>% filter(!is.na(spatial_id)) # those with spatial counterparts

  # scale data on all features
  multiome <- ScaleData(multiome, features = rownames(multiome))
  spatial <- ScaleData(spatial, features = rownames(spatial))

  shared_gene_SumExp <- data.frame(SumExp = colSums(multiome@assays[[multiome_assay]]@scale.data[shared_genes, ])) %>% rownames_to_column(var = 'multiome_id')

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.spatial.coord.summed.expression.pdf'), width = 5, height = 4)}
  print(outs.metadata %>% dplyr::select(multiome_id, x, y) %>% left_join(shared_gene_SumExp, by = 'multiome_id') %>%
          ggplot( aes(y, x, color = SumExp)) + geom_point(size = 0.8) + theme_classic() + scale_color_viridis_c( direction = -1) +
          ggtitle('Summed expression of shared genes in multiome') + theme(plot.title = element_text(hjust = 0.5)))
  if(!is.null(plot_dir)){dev.off()}


  # for each cell - each point represents a cell
  cor_pr_cell <- list(); cor_sp_cell <- list()
  for(i in 1:nrow(ids)){
    multiome_exp <- as.matrix(multiome@assays[[multiome_assay]]@scale.data[shared_genes, ids$multiome_id[i] ])
    spatial_exp <- as.matrix(spatial@assays[[spatial_assay]]@scale.data[shared_genes,ids$spatial_id[i] ])
    cor_pr_cell[i] <- cor(multiome_exp, spatial_exp, method = 'pearson', use = "complete.obs")
    cor_sp_cell[i] <- cor(multiome_exp, spatial_exp, method = 'spearman', use = "complete.obs")
  }

  cor_cell <- cbind(t(data.frame(cor_pr_cell)), t(data.frame(cor_sp_cell)))
  colnames(cor_cell) <- c('pearson', 'spearman')

  # show correlation on spatial slide
  cor_cell_ids <- cbind(ids, cor_cell)
  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.scatterplot.summed.expression.wrt.spearman.correlation.pdf'), width = 5, height = 4)}
  print(cor_cell_ids %>% left_join(shared_gene_SumExp, by = 'multiome_id') %>%
          ggplot( aes(SumExp, spearman)) + geom_point(size = 0.5) + theme_classic() +
          labs(x = 'Summed Expression of\nMultiome-Spatial Shared Genes', y = 'Spearman correlation\nbetween Multiome-Spatial') +
          ggtitle('Correlate summed expression with projection accuracy') + theme(plot.title = element_text(hjust = 0.5)))
  if(!is.null(plot_dir)){dev.off()}

  outs <- outs.metadata %>% dplyr::select(multiome_id, x, y) %>% left_join(shared_gene_SumExp, by = 'multiome_id')
  return(outs)


}


#' @title Plot nearest neighbor distance
#' @description Diagnostic plot showing distance of closest spatial counterpart with multiome cells
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return Diagnostic plots describing how close spatial counterparts were to each multiome cells
#' @export
PlotNNDistance <- function(outs.metadata, plot_dir = NULL){

  # diagnostic plots
  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Diagnostic.scatterplot.nn.k.vs.distance.pdf'), width = 4, height = 3.5)}
  print(outs.metadata %>%
          ggplot(aes(x=as.numeric(nn_k), y=nn_dist)) +
          geom_point() +
          stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
          theme_classic() +  xlab("Closest spatial neighbor's K value") + ylab('Distance with Multiome') )
  if(!is.null(plot_dir)){dev.off()}

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Diagnostic.spatial.coord.nn.k.pdf'), width = 5, height = 4)}
  print(outs.metadata %>% mutate(nn_k = as.numeric(nn_k)) %>%
          ggplot( aes(y, x, color = nn_k)) + geom_point(size = 0.8)  +
          theme_classic() + scale_color_distiller(palette = "YlGnBu", direction = 1) +
          labs(color = "Closest spatial\nneighbor's K")
  )
  if(!is.null(plot_dir)){dev.off()}


  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Diagnostic.spatial.coord.nn.distance.pdf'), width = 5, height = 4)}
  print(outs.metadata %>%
          ggplot( aes(y, x, color = nn_dist)) + geom_point(size = 0.8)  +
          theme_classic() +scale_color_distiller(palette = "YlGnBu", direction = 1) +
          labs(color = "Distance\nwith Multiome")
  )
  if(!is.null(plot_dir)){dev.off()}


}

#' @title Plot spatial assignment
#' @description Diagnostic plot showing number of multiome cells assigned to each spatial coordinate
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return A dataframe showing number of multiome cells assigned for unique spatial coordinate.
#' Plots histogram and scatterplot describing the number of multiome cells assigned to each spatial coordinate
#' @export
PlotSpatialAssignment <- function(outs.metadata, plot_dir = NULL){
  outs.metadata$coord <- paste0(round(outs.metadata$x, 2), '_', round(outs.metadata$y, 2))

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Diagnostic.histogram.multiome.assignment.pdf'), width = 4, height = 3.5)}
  print(as.data.frame(table(outs.metadata$coord)) %>%
          ggplot(aes(x=Freq)) +
          geom_histogram(color="white", fill="#7f53a6") + theme_classic() +
          labs(x = '# of Multiome cells assigned per spatial coordinate', y = '# of Spatial coordinates') )
  if(!is.null(plot_dir)){dev.off()}

  freq <- outs.metadata %>% left_join(as.data.frame(table(outs.metadata$coord)) , by = c('coord' = 'Var1'))

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Diagnostic.spatial.coord.multiome.assignment.pdf'), width = 5, height = 4)}
  print(freq %>%
          ggplot( aes(y, x, color = Freq)) + geom_point(size = 0.8)  +
          theme_classic() + scale_color_distiller(palette = "PuRd", direction = 1) +
          #scale_color_viridis_c(option = "inferno", direction = -1) +
          labs(color = "# of Multiome\ncells located")
  )
  if(!is.null(plot_dir)){dev.off()}


  return(freq)
}
