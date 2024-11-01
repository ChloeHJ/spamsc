#' @title PlotCorrelation
#'
#' @description A plotting function showing correlation of gene expression from spatial and multiome data
#'
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param spatial_assay assay within spatial object shared with multiome data to integrate
#' @param multiome_assay assay within multiome object shared with spatial data to integrate
#' @param log10 whether to plot correlation on after log10 transformation, default = TRUE
#' @param plot_dir Directory to save plots. If NULL, print plots.
#'
#' @return A scatterplot showing summed gene counts from spatial and multome data
#' @export
PlotCorrelation <- function(spatial, multiome, spatial_assay = 'Xenium', multiome_assay = 'RNA', log10 = TRUE, plot_dir = NULL){

  # correlate multiome & spatial expression
  DefaultAssay(spatial) <- spatial_assay
  DefaultAssay(multiome) <- multiome_assay

  spatial <- CreateSeuratObject(spatial[[spatial_assay]], assay =spatial_assay)
  multiome <- CreateSeuratObject(multiome[[multiome_assay]], assay = multiome_assay)

  shared_genes <- intersect(rownames(spatial), rownames(multiome))
  print(paste0('There are ', length(shared_genes), ' shared genes between spatial and multiome data for integration'))
  spatial <- spatial[shared_genes, ]
  multiome <- multiome[shared_genes, ]

  # correlate expression (scatterplot by counts)
  spatial_gene_counts <- rowSums(spatial@assays[[spatial_assay]]@data)
  multiome_gene_counts <- rowSums(multiome@assays[[multiome_assay]]@data)

  sp_gene_count <- data.frame(spatial_gene_counts) %>% rownames_to_column(var = 'gene')
  sc_gene_count <- data.frame(multiome_gene_counts) %>% rownames_to_column(var = 'gene')

  x_lab <- paste0('Gene counts in Multiome ', multiome_assay)
  y_lab <- paste0('Gene counts in Spatial ', spatial_assay)

  gene_count_dt <- sp_gene_count %>% inner_join(sc_gene_count, by = 'gene')
  if(log10 == TRUE){
    if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/0_scatterplot.correlate.multiome.spatial.expression.pdf'), width = 5, height = 4)}
    print (gene_count_dt %>%
             ggplot(aes(x=multiome_gene_counts, y=spatial_gene_counts)) +
             geom_point() +
             geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
             theme_classic() + scale_x_continuous(trans='log10') +
             scale_y_continuous(trans='log10') + xlab(x_lab) + ylab(y_lab) + labs(title = paste0(length(shared_genes),' Shared genes')))
    if(!is.null(plot_dir)){dev.off()}

  }else{
    if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/0_scatterplot.correlate.multiome.spatial.expression.pdf'), width = 5, height = 4)}
    print (gene_count_dt %>%
             ggplot(aes(x=multiome_gene_counts, y=spatial_gene_counts)) +
             geom_point() +
             geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
             theme_classic() + xlab(x_lab) + ylab(y_lab)  + labs(title = paste0(length(shared_genes),' Shared genes')))
    if(!is.null(plot_dir)){dev.off()}

  }


  return(gene_count_dt)

}


#' @title MergeDatasets
#' @description Integrate spatial and multiome datasets
#'
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param spatial_assay assay within spatial object shared with multiome data to integrate
#' @param multiome_assay assay within multiome object shared with spatial data to integrate
#' @param spatial_cluster column name in meta.data of spatial Seurat object containing cluster annotation of spatial data
#' @param multiome_cluster column name in meta.data of multiome Seurat object containing cluster annotation of multiome data
#' @param coord cell x 2 matrix containing spatial coordinate of cells
#' @param lambda Ridge regression penalty parameter for Harmony integration.
#' Default lambda=NULL for automatic lambda estimation.
#' Lambda must be strictly positive. Smaller values result in more aggressive correction.
#' @param gene_counts A threshold to remove genes based on summed gene counts. Default: 0.
#' i.e., only retain genes having summed gene counts > gene_counts.
#' @param cell_counts A threshold to remove cells based on summed cell counts. Default: 0.
#' i.e., only retain cells having summed cell counts > cell_counts
#' @param plot_dir Directory to save plots. If NULL, print plots.
#'
#' @return A list of integrated Seurat object from spatial and multiome data
#' and coordinate after filtering for cells with low counts
#' @export
MergeDatasets <- function (spatial, multiome, spatial_assay = "Xenium", multiome_assay = "RNA",
          spatial_cluster = "RNACluster", multiome_cluster = "RNA_Label",
          coord, lambda = NULL, gene_counts = 0, cell_counts = 0, plot_dir = NULL)
{

  DefaultAssay(spatial) <- spatial_assay
  DefaultAssay(multiome) <- multiome_assay

  # pre-processing
  multiome <- multiome[which(rowSums(multiome@assays[[multiome_assay]]@counts) >  gene_counts), ]
  spatial <- spatial[which(rowSums(spatial@assays[[spatial_assay]]@counts) > gene_counts), ]
  shared_genes <- intersect(rownames(spatial), rownames(multiome))
  print(paste0("There are ", length(shared_genes), " shared genes between spatial and multiome data for integration"))
  spatial_data <- spatial[shared_genes, ]
  multiome_data <- multiome[shared_genes, ]

  multiome_data <- multiome_data[, which(colSums(multiome_data@assays[[multiome_assay]]@counts) >   cell_counts)]
  spatial_data <- spatial_data[, which(colSums(spatial_data@assays[[spatial_assay]]@counts) >  cell_counts)]
  coord <- coord[which(colSums(spatial_data@assays[[spatial_assay]]@counts) > cell_counts), ]

  # add cluster labels and datasets
  multiome_data$PDataset <- "Multiome"
  multiome_data$PCluster <- multiome_data[[multiome_cluster]]
  spatial_data$PDataset <- "Spatial"
  spatial_data$PCluster <- spatial_data[[spatial_cluster]]

  # just in case seurats were from previous versions, to make compatible with seurat v5 worfklow
  multiome_assay5 <- as(multiome_data[[multiome_assay]], Class = "Assay5")
  multiome_data <- CreateSeuratObject(multiome_assay5, meta.data = multiome_data@meta.data, assay = 'RNA')
  spatial_assay5 <- as(spatial_data[[spatial_assay]], Class = "Assay5")
  spatial_data <- CreateSeuratObject(spatial_assay5, meta.data = spatial_data@meta.data, assay = 'RNA')

  # merge by seurat v5 workflow
  merge <- merge(spatial_data, multiome_data)
  merge@active.assay <- "RNA"
  merge <- JoinLayers(merge, assay = "RNA")

  # batch correct
  merge[["RNA"]] <- split(merge[["RNA"]], f = merge$PDataset)
  merge <- SCTransform(merge, assay = "RNA")
  merge <- RunPCA(object = merge)
  merge <- IntegrateLayers(merge, method = HarmonyIntegration,
                           orig = "pca", new.reduction = "harmony", dims = 1:30,
                           features = rownames(merge), normalization.method = "SCT", lambda = lambda,
                           verbose = T)
  merge <- RunUMAP(merge, dims = 1:10, reduction = "harmony",  reduction.name = "umap.harmony")
  if (!is.null(plot_dir)) {
    cairo_pdf(filename = paste0(plot_dir, "/1_umap.integrated.data.color.by.dataset.pdf"),
              width = 6, height = 5)
  }
  print(DimPlot(merge, group.by = "PDataset", raster = FALSE,
                reduction = "umap.harmony") + ggtitle("Integrated data, color by Dataset") +
          theme(plot.title = element_text(hjust = 0.5)))
  if (!is.null(plot_dir)) {
    dev.off()
  }
  if (!is.null(plot_dir)) {
    cairo_pdf(filename = paste0(plot_dir, "/1_umap.integrated.data.color.by.cluster.pdf"),
              width = 9, height = 5)
  }
  print(DimPlot(merge, group.by = "PCluster", raster = FALSE,
                reduction = "umap.harmony") + ggtitle("Integrated data, color by Cluster") +
          theme(plot.title = element_text(hjust = 0.5)))
  if (!is.null(plot_dir)) {
    dev.off()
  }
  if (!is.null(plot_dir)) {
    cairo_pdf(filename = paste0(plot_dir, "/1_umap.integrated.data.color.by.spatial.cluster.pdf"),
              width = 7.5, height = 5)
  }
  print(DimPlot(merge, group.by = spatial_cluster, raster = FALSE,
                reduction = "umap.harmony") + ggtitle("Integrated data, color by spatial clusters") +
          theme(plot.title = element_text(hjust = 0.5)))
  if (!is.null(plot_dir)) {
    dev.off()
  }

  coord_dt <-  data.frame(spatial_id = colnames(spatial_data),coord)
  outs <- list(merge = merge, coord_dt = coord_dt)

  return(outs)
}

#' @title RunProjection
#' @description Core function to integrate spatial and multiome data, and assign spatial coordinate to multiome cells
#'
#' @param merge A merged Seurat object, output from `MergeDatasets()` function, merging spatial and multiome datasets
#' @param coord_dt cell x 3 data frame containing cell names and spatial coordinate of cells, output from `MergeDatasets()`
#' @param knn_k a cutoff for the top k k-nearest neighbors to find a spatial counterpart for. Default: 100.
#' i.e., if there is no spatial counterpart in the top 100 neighbors, no spatial coordinate
#' is assigned to the respective multiome cell.
#' coordinate is assigned for the respective multiome cell.
#' @param n_harmony_dim Number of harmony dimension to use to compute nearest neighbour. Default: 35
#' @param multiome_cluster column name in meta.data of merge Seurat object to visualise after projection
#' @param plot_dir Directory to save plots. If NULL, print plots.
#'
#' @return  a metadata containing spatial coordinate of multliome cells
#' @export
RunProjection <- function(merge, coord_dt, knn_k = 100, n_harmony_dim = 35,  multiome_cluster = "RNA_Label",  plot_dir = NULL){

  # distance on harmony umap
  print('Computing nearest neighbour...')
  har_emb <- merge@reductions$harmony@cell.embeddings[, 1:n_harmony_dim] #knn based on harmony umap embedding
  knn_k <- knn_k #[P]k threshold to get the best counterpart ie. if the best 100 neighbours are all multiome, then no spatial ids are retrieved
  nn <- kNN(har_emb, k = knn_k)

  i <- 10 #NOT PARAMETER arbitraty value for visualisation
  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/2_umap.nearest.neighbours.harmony.pdf'), width = 6, height = 5)}
  print(DimPlot(object = merge, cells.highlight =list(cell = colnames(merge)[i], neighbours = c(colnames(merge)[nn$id[i,]]) ),
                cols.highlight = c("orange", 'red'), cols = "gray", raster = FALSE, reduction = 'harmony') +
          ggtitle('Nearest top k cells, in harmony')+ theme(plot.title = element_text(hjust = 0.5)) )
  if(!is.null(plot_dir)){dev.off()}

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/2_umap.nearest.neighbours.umap.harmony.pdf'), width = 6, height = 5)}
  print(DimPlot(object = merge, cells.highlight =list(cell = colnames(merge)[i], neighbours = c(colnames(merge)[nn$id[i,]]) ),
                cols.highlight = c("orange", 'red'), cols = "gray", raster = FALSE, reduction = 'umap.harmony') +
          ggtitle('Nearest top k cells, in umap.harmony')+ theme(plot.title = element_text(hjust = 0.5)) )
  if(!is.null(plot_dir)){dev.off()}


  # get best spatial counterpart of multimodal
  spaital_ids <- which(merge$PDataset == 'Spatial')
  multiome_ids <- which(merge$PDataset == 'Multiome')
  nn_id_mtx <- nn$id[multiome_ids, ] # get multiome only
  nn_dist_mtx <- nn$dist[multiome_ids, ]

  #get best spatial counterpart indices
  print('Finding the best spatial counterpart...')
  best_spatial_inds <- c()
  for(i in 1:nrow(nn_id_mtx)){
    row <- nn_id_mtx[i, ]
    row[which(!row %in% spaital_ids)] <- NA
    nn_id_mtx[i, ] <- row
    best_spatial_inds <- c(best_spatial_inds, row[which(!is.na(row))[1]])
  }

  # output dataframe
  sp_counterparts <- data.frame(
    multiome_id = rownames(nn_id_mtx),
    nn_k = names(best_spatial_inds),
    spatial_ind = best_spatial_inds
    ) %>%
    mutate(spatial_id = rownames(nn$id)[spatial_ind] )

  # get distance
  for(row in 1:nrow(sp_counterparts)){
    k <- as.numeric(sp_counterparts$nn_k[row])
    sp_counterparts$nn_dist[row] <- nn_dist_mtx[i, k]
  }

  # get location of spatial
  sp_counterparts <- sp_counterparts %>% left_join(coord_dt, by = c('spatial_id'))

  # add to multiome metadata
  outs.metadata <- multiome@meta.data %>% rownames_to_column(var =  'multiome_id') %>%
    left_join(sp_counterparts, by = 'multiome_id')

  if(!is.null(plot_dir)){cairo_pdf(filename = paste0(plot_dir, '/3_spatial.location.of.multiome.pdf'), width = 7, height = 5)}
  print( outs.metadata %>%  ggplot( aes(y, x, color = outs.metadata[[multiome_cluster]])) + geom_point(size = 0.7) + theme_classic() +
           ggtitle('Spatial location of Multiome cells, color by Multiome clusters') + theme(plot.title = element_text(hjust = 0.5)) +
           guides(color=guide_legend(title="Multiome cluster")))
  if(!is.null(plot_dir)){dev.off()}

  print('Done!')
  return(outs.metadata)
}



