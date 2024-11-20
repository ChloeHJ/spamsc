#' @title VisualisePairedModalityOneCounterpart
#' @description Visualize gene expression from twp multiome assays on spatial coordinates
#' when there is only one counterpart between assay 1 (e.g., RNA assay) and assay 2 (e.g., protein assay).
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param lookup_table A table containing features of the multiome assay used for integration (named 'gene')
#' and the corresponding features in the assay to project (named 'paired'). If the feature names are the same,
#' duplicate the column with the names 'gene' and 'paired'.
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param multiome_assay assay within multiome object shared with spatial data
#' @param multiome_project_assay assay within multiome object to project
#' @param gene gene of interest
#' @param pt.size size of points on plot. Default = 0.5
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return A scatterplot showing gene expression in the shared and projected assays on spatial coordinates.
#' @export
VisualisePairedModalityOneCounterpart <- function(outs.metadata, lookup_table, multiome, spatial,
                                                  gene,  multiome_assay = 'RNA', multiome_project_assay = 'ADT',
                                                  pt.size = 0.5, plot_dir = NULL){

  if(gene %in% rownames(spatial)){
    print(ImageFeaturePlot(spatial, features = gene, size = 1.2) & scale_fill_viridis_c(direction = 1) )
  }else{print('Gene not found on spatial slide')}

  genes_exp <- tryCatch({
    as.matrix(multiome@assays[[multiome_assay]]@data[gene, ]) %>%
      as.data.frame() %>% rownames_to_column(var = 'multiome_id')
  }, error = function(e) {
    as.matrix(multiome@assays[[multiome_assay]]@counts[gene, ]) %>%
      as.data.frame() %>% rownames_to_column(var = 'multiome_id')
  })

  colnames(genes_exp) <- c('multiome_id', 'V1')
  genes_exp$V1 <- scale(genes_exp$V1)

  p1 <- outs.metadata %>% dplyr::select(multiome_id, x, y) %>% left_join(genes_exp, by = 'multiome_id') %>%
    ggplot( aes(y, x, color = V1)) + geom_point(size = pt.size) + theme_classic() +
    scale_color_viridis(direction = -1)  + labs(color = gene) + ggtitle(paste0(multiome_assay, ' : ', gene)) +
    labs(color='Scaled\nExpression')

  # paired modality
  paired_gene <- lookup_table$paired[which(lookup_table$gene == gene)]
  if(identical(paired_gene, character(0))){
    print('Gene not found in lookup table')
    print(p1)
  }else{

    genes_exp <- tryCatch({
      as.matrix(multiome@assays[[multiome_project_assay]]@data[paired_gene, ]) %>%
        as.data.frame() %>% rownames_to_column(var = 'multiome_id')
    }, error = function(e) {
      as.matrix(multiome@assays[[multiome_project_assay]]@counts[paired_gene, ]) %>%
        as.data.frame() %>% rownames_to_column(var = 'multiome_id')
    })


    colnames(genes_exp) <- c('multiome_id', 'V1')
    genes_exp$V1 <- scale(genes_exp$V1)

    p2 <- outs.metadata %>% dplyr::select(multiome_id, x, y) %>%
      left_join(genes_exp, by = 'multiome_id') %>%
      ggplot( aes(y, x, color = V1)) + geom_point(size = pt.size) + theme_classic() +
      scale_color_viridis(direction = -1) + labs(color = paired_gene) +
      ggtitle(paste0(multiome_project_assay, ' : ', gene)) + labs(color='Scaled\nExpression')

    if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Visualise.spatial.coord.', gene, '.', multiome_project_assay, '.pdf'), width = 7, height = 2.5)}
    print(p1 + p2)
    if(!is.null(plot_dir)){dev.off()}

  }

  return(p1 + p2)

}

#' @title VisualisePairedModality
#' @description Visualize gene expression from two multiome assays on spatial coordinates
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param lookup_table A table containing features of the multiome assay used for integration (named 'gene')
#' and the corresponding features in the assay to project (named 'paired'). If the feature names are the same,
#' duplicate the column with the names 'gene' and 'paired'.
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param multiome_assay assay within multiome object shared with spatial data
#' @param multiome_project_assay assay within multiome object to project
#' @param gene gene of interest
#' @param summary_mode When there are multiple counterparts between assay 1 (e.g., RNA assay)
#' and assay 2 (e.g., ATAC peaks), there are modes to summarize the expression. `all` plots all
#' features corresponding to the gene, while `mean` and `median` summarize the expression into a single plot.
#' @param pt.size size of points on plot. Default = 0.5
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return A scatterplot showing gene expression in the shared and projected assays on spatial coordinates.
#' @export
VisualisePairedModality <- function(outs.metadata, lookup_table, multiome, spatial, gene,
                                    multiome_assay = 'RNA', multiome_project_assay = 'ADT',
                                    summary_mode = c('all', 'mean', 'median'), pt.size = 0.5, plot_dir = NULL){

  paired_gene <- lookup_table$paired[which(lookup_table$gene == gene)]
  if(length(paired_gene) == 0){
    'Gene not found in projected modality'
  }else if(length(paired_gene)  == 1){
    print('One counterpart found in paired modality')
    VisualisePairedModalityOneCounterpart(outs.metadata = outs.metadata,
                                          lookup_table = lookup_table,
                                          multiome = multiome,
                                          spatial = spatial,
                                          gene =gene,
                                          multiome_assay = multiome_assay,
                                          multiome_project_assay = multiome_project_assay,
                                          pt.size = pt.size, plot_dir = plot_dir)

  }else{

    print('More than one counterpart found in paired modality')

    # spatial
    if(gene %in% rownames(spatial)){
      print(ImageFeaturePlot(spatial, features = gene, size = 1.2) & scale_fill_viridis_c(direction = 1) )
    }else{print('Gene not found on spatial slide')}


    genes_exp <- tryCatch({
      as.matrix(multiome@assays[[multiome_assay]]@data[gene, ]) %>%
        as.data.frame() %>% rownames_to_column(var = 'multiome_id')
    }, error = function(e) {
      as.matrix(multiome@assays[[multiome_assay]]@counts[gene, ]) %>%
        as.data.frame() %>% rownames_to_column(var = 'multiome_id')
    })

    colnames(genes_exp) <- c('multiome_id', 'V1')
    genes_exp$V1 <- scale(genes_exp$V1)
    p1 <- outs.metadata %>% dplyr::select(multiome_id, x, y) %>% left_join(genes_exp, by = 'multiome_id') %>%
      ggplot( aes(y, x, color = V1)) + geom_point(size = pt.size) + theme_classic() +
      scale_color_viridis(direction = -1)  + labs(color = gene) +  ggtitle(paste0(multiome_assay, ' : ', gene)) +
      labs(color='Scaled\nExpression')

    # paired modality
    if(summary_mode == 'all'){

      genes_exp <- tryCatch({
        as.matrix(multiome@assays[[multiome_project_assay]]@data[paired_gene, ]) %>%
          as.data.frame() %>% rownames_to_column(var = 'multiome_id')
      }, error = function(e) {
        as.matrix(multiome@assays[[multiome_project_assay]]@counts[paired_gene, ]) %>%
          as.data.frame() %>% rownames_to_column(var = 'multiome_id')
      })

      all_gene_exp <- melt(genes_exp, id.var = 'multiome_id')
      colnames(all_gene_exp) <- c('multiome_id', 'variable' ,'value')
      all_gene_exp <- all_gene_exp %>% group_by(variable) %>% mutate(value = scale(value))
      p2 <- outs.metadata %>% dplyr::select(multiome_id, x, y) %>%
        left_join(all_gene_exp, by = 'multiome_id') %>%
        ggplot( aes(y, x, color = value)) + geom_point(size = pt.size) + theme_classic() +
        scale_color_viridis(direction = -1) +  facet_wrap(~variable) +  ggtitle(paste0(multiome_project_assay, ' : ', gene)) +
        labs(color='Scaled\nExpression')

      if(!is.null(plot_dir)){ cairo_pdf(filename = paste0(plot_dir, '/Visualise.spatial.coord.all.', gene, '.', multiome_project_assay, '.pdf'), width = 7, height = 2.5)}
      print(p1 + p2)
      if(!is.null(plot_dir)){dev.off()}


    }

    if(summary_mode == 'mean'){

      genes_exp <- tryCatch({
        as.matrix(multiome@assays[[multiome_project_assay]]@data[paired_gene, ]) %>%
          as.data.frame() %>% rownames_to_column(var = 'multiome_id')
      }, error = function(e) {
        as.matrix(multiome@assays[[multiome_project_assay]]@counts[paired_gene, ]) %>%
          as.data.frame() %>% rownames_to_column(var = 'multiome_id')
      })

      genes_exp <- as.data.frame(cbind(genes_exp[, 1], rowMeans(genes_exp[, 2:ncol(genes_exp)])))
      colnames(genes_exp) <- c('multiome_id', 'V1')
      genes_exp$V1 <- scale(genes_exp$V1)

      p2 <- outs.metadata %>% dplyr::select(multiome_id, x, y) %>% left_join(genes_exp, by = 'multiome_id') %>%
        mutate(V1 = as.numeric(V1)) %>%
        ggplot( aes(y, x, color = V1)) + geom_point(size = pt.size) + theme_classic() +
        scale_color_viridis(direction = -1) + labs(color = paste0(gene, '-Mean')) +
        ggtitle(paste0(multiome_project_assay, ' : ', gene)) +
        labs(color='Scaled\nExpression')

      if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Visualise.spatial.coord.mean.', gene, '.', multiome_project_assay, '.pdf'), width = 7, height = 2.5)}
      print(p1 + p2)
      if(!is.null(plot_dir)){dev.off()}

    }

    if(summary_mode == 'median'){

      genes_exp <- tryCatch({
        as.matrix(multiome@assays[[multiome_project_assay]]@data[paired_gene, ]) %>%
          as.data.frame() %>% rownames_to_column(var = 'multiome_id')
      }, error = function(e) {
        as.matrix(multiome@assays[[multiome_project_assay]]@counts[paired_gene, ]) %>%
          as.data.frame() %>% rownames_to_column(var = 'multiome_id')
      })

      genes_exp <- as.data.frame(cbind(genes_exp[, 1], rowMedians(as.matrix(genes_exp[, 2:ncol(genes_exp)]))))
      colnames(genes_exp) <- c('multiome_id', 'V1')
      genes_exp$V1 <- scale(genes_exp$V1)
      p2 <- outs.metadata %>% dplyr::select(multiome_id, x, y) %>%
        left_join(genes_exp, by = 'multiome_id') %>% mutate(V1 = as.numeric(V1)) %>%
        ggplot( aes(y, x, color = V1)) + geom_point(size = pt.size) + theme_classic() +
        scale_color_viridis(direction = -1) +
        labs(color = paste0(gene, '-Median')) +  ggtitle(paste0(multiome_project_assay, ' : ', gene)) +
        labs(color='Scaled\nExpression')

      if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Visualise.spatial.coord.median.', gene, '.', multiome_project_assay, '.pdf'), width = 7, height = 2.5)}
      print(p1 + p2)
      if(!is.null(plot_dir)){dev.off()}


    }

  }

  return(p1 + p2)

}


#' @title CompareClusters
#' @description Visualize the proportion of cells in multiome and spatial clusters.
#' Visualize the proportion of cells in multiome clusters that lie in spatial clusters, and vice versa.
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param spatial_cluster column name in meta.data of spatial Seurat object containing cluster annotation of spatial data
#' @param multiome_cluster column name in meta.data of multiome Seurat object containing cluster annotation of multiome data
#' @param width width of the plot to save
#' @param height height of the plot to save
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#' @return A heatmap showing the absolute and relative overlap of cells between two spatial and multiome clusters
#' @export
CompareClusters <- function(outs.metadata, multiome, spatial,
                            spatial_cluster = 'RNACluster', multiome_cluster = 'RNA_Clusters_HC',
                            width = 4, height = 3, plot_dir = NULL){

  id_match <- outs.metadata %>% select(multiome_id, spatial_id)
  length(which(is.na(id_match$spatial_id))) #multiome cells don't have spatial counterparts
  ids <- id_match %>% filter(!is.na(spatial_id)) # those with spatial counterparts

  clusters <- ids %>%
    left_join(data.frame(Multiome_Label = multiome[[multiome_cluster]] ) %>% rownames_to_column(var = 'multiome_id'), by = 'multiome_id') %>%
    left_join(data.frame(Spatial_Label = spatial[[spatial_cluster]]) %>% rownames_to_column(var = 'spatial_id'), by ='spatial_id' )
  colnames(clusters) <- c('multiome_id', 'spatial_id', 'Multiome_Label', 'Spatial_Label')


  table <- table(clusters$Multiome_Label, clusters$Spatial_Label)
  table_lg <- table %>% melt() %>% mutate(Var2 = as.character(Var2))
  colnames(table_lg) <- c('Multiome_Label', 'Spatial_Label', 'freq')

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.compare.multiome.spatial.clusters.freq.pdf'), width = width, height = height)}
  print(ggplot(table_lg, aes(Multiome_Label, Spatial_Label, fill= freq)) +
          geom_tile() + scale_fill_gradient(low="white", high="red") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) )
  if(!is.null(plot_dir)){dev.off()}

  p <- table_lg %>% group_by(Spatial_Label) %>% mutate(sum_spatial = sum(freq)) %>%
    mutate(`% Spatial` = freq/sum_spatial) %>%
    ggplot(aes(Multiome_Label, Spatial_Label, fill= `% Spatial`)) +
    geom_tile() + scale_fill_gradient(low="white", high="red") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.compare.multiome.spatial.clusters.pct.spatial.pdf'), width = width, height = height)}
  print( p )
  if(!is.null(plot_dir)){dev.off()}

  p <- table_lg %>% group_by(Multiome_Label) %>% mutate(sum_multiome = sum(freq)) %>%
    mutate(`% Multiome` = freq/sum_multiome) %>%
    ggplot(aes(Multiome_Label, Spatial_Label, fill= `% Multiome`)) +
    geom_tile() + scale_fill_gradient(low="white", high="red") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/Evaluate.compare.multiome.spatial.clusters.pct.multiome.pdf'), width = width, height = height)}
  print( p )
  if(!is.null(plot_dir)){dev.off()}

  return(p)

}


#' @title VisualiseSingleModality
#' @description Visualize gene expression from any multiome assays on spatial coordinates
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param multiome A multiome Seurat object
#' @param multiome_plot_assay assay within multiome object to plot
#' @param gene gene of interest
#' @param pt.size size of points on plot. Default = 0.5
#' @param order order cells by the gene expression
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return A scatterplot showing gene expression from any mulitome assay on spatial coordinates.
#' @export
VisualiseSingleModality <- function(outs.metadata, multiome, gene = NA,
                                    multiome_plot_assay = 'GeneScore', pt.size = 0.4, order = FALSE, plot_dir = NULL,
                                    plot_height = 3, plot_width = 4){

  # original modality
  DefaultAssay(multiome) <- multiome_plot_assay;
  outs.metadata$multiome_id <- rownames(outs.metadata)

  genes_exp <- tryCatch({
    as.matrix(multiome@assays[[multiome_plot_assay]]@counts[gene, ]) %>%
      as.data.frame() %>% rownames_to_column(var = 'multiome_id')
  }, error = function(e) {
    as.matrix(multiome@assays[[multiome_plot_assay]]@data[gene, ]) %>%
      as.data.frame() %>% rownames_to_column(var = 'multiome_id')
  })
  genes_exp$V1 <- scale(genes_exp$V1)
  colnames(genes_exp) <- c('multiome_id', 'V1')
  data <- outs.metadata %>% dplyr::select(multiome_id, x, y) %>% left_join(genes_exp, by = 'multiome_id')
  if(order == TRUE){data <- data %>% arrange(V1)}

  p <- data %>%
    ggplot( aes(y, x, color = V1)) + geom_point(size = pt.size) + theme_classic() +
    scale_color_distiller(palette = "Spectral") + labs(color = gene)  +
    ggtitle(paste0(multiome_plot_assay, ': ', gene)) + labs(color='Scaled\nExpression')
  if(!is.null(plot_dir)){
    if(order == TRUE){cairo_pdf(filename = paste0(plot_dir, '/Visualise.spatial.coord.', gene, '.', multiome_plot_assay, '.ordered.pdf'), width = plot_width, height = plot_height)}else{
      cairo_pdf(filename = paste0(plot_dir, '/Visualise.spatial.coord.', gene, '.', multiome_plot_assay, '.pdf'), width = plot_width, height = plot_height)
    }
      }
  print( p)
  if(!is.null(plot_dir)){dev.off()}

  return(p)

}


#' @title VisualiseFeatures
#' @description Visualize feature expression from meta.data on spatial coordinates
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param multiome A multiome Seurat object
#' @param feature feature of interest e.g. RNA counts, cluster labels
#' @param pt.size size of points on plot. Default = 0.5
#' @param order order cells by the gene expression
#' @param plot_dir Directory to save plots. If NULL, print plots. Example: '/directory' (instead of '/directory/')
#'
#' @return A scatterplot showing feature on spatial coordinates.
#' @export
VisualiseFeatures <- function(outs.metadata, multiome, feature = NA, pt.size = 0.4, order = FALSE, plot_dir = NULL){


  genes_exp <- as.matrix(multiome[[feature]]) %>% as.data.frame() %>% rownames_to_column(var = 'multiome_id')
  colnames(genes_exp) <- c('multiome_id', 'V1')
  data <- outs.metadata %>% dplyr::select(multiome_id, x, y) %>%  left_join(genes_exp, by = 'multiome_id')
  if(order == TRUE){data <- data %>% arrange(V1)}

  if(class(genes_exp$V1) == 'numeric'){
    p <- data %>%
      ggplot( aes(y, x, color = V1)) + geom_point(size = pt.size) + theme_classic() +
      labs(color = feature) + scale_color_viridis(direction = -1)
  }else{
    p <- data %>%
      ggplot( aes(y, x, color = V1)) + geom_point(size = pt.size) +
      theme_classic() + labs(color = feature)
  }

  if(!is.null(plot_dir)){
    if(order == TRUE){ cairo_pdf(filename = paste0(plot_dir, '/Visualise.spatial.coord.', feature, '.ordered.pdf'), width = 3.5, height = 3)}else{
      cairo_pdf(filename = paste0(plot_dir, '/Visualise.spatial.coord.', feature, '.pdf'), width = 5, height = 4)}
    }
  print(p)
  if(!is.null(plot_dir)){dev.off()}

  return(p)
}
