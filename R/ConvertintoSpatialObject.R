#' @title Convert Multiome to SpatialObject
#' @description Processing function to convert a multiome Seurat object into a spatial Seurat object by
#' adding spatial coordinates and image object
#'
#' @param outs.metadata A output from `run_projection()`, containing spatial coordinate of multiome cells
#' @param spatial A spatial Seurat object
#' @param multiome A multiome Seurat object
#' @param spatial_assay assay within spatial object shared with multiome data to integrate
#' @param image am image (fov) from spatial data to include in the multiome spatial object
#'
#' @return A spatial Seurat object
#' @export
MultiometoSpatialObject <- function(outs.metadata, spatial, multiome, spatial_assay = 'MERSCOPE', image = 'AAH_COLON_6_MONTHS'){

  metadata <- outs.metadata %>% select(multiome_id, spatial_id)

  # subset spatial to matched cells only
  spatial$spatial_id <- colnames(spatial)
  spatial.subset <- subset(spatial, subset = spatial_id %in% metadata$spatial_id)
  spatial.subset@meta.data <- spatial.subset@meta.data  %>%
    rownames_to_column(var = 'rowname') %>%
    left_join(metadata, by = 'spatial_id') %>%
    group_by(spatial_id) %>% slice_head(n = 1) %>% ungroup() %>% # currently only subsetting first 1-1 match
    column_to_rownames(var = 'rowname')
  spatial.subset <- RenameCells(spatial.subset, new.names = spatial.subset$multiome_id)

  # subset multiome to matched cells only
  multiome$multiome_id <- colnames(multiome)
  multiome.subset <- subset(multiome, subset = multiome_id %in% colnames(spatial.subset))

  # add spatial assay and image into multiome data
  multiome.subset[['Spatial']]  <- CreateAssayObject(counts = spatial.subset@assays[[spatial_assay]]@counts[, match(colnames(multiome.subset), colnames(spatial.subset))],
                                                     assay = 'Spatial')
  multiome.subset@images[[image]] <- spatial.subset@images[[image]][, match(colnames(multiome.subset), colnames(spatial.subset))]
  multiome.subset@meta.data <- multiome.subset@meta.data %>%
    rownames_to_column(var = 'rowname') %>%
    left_join(outs.metadata %>% select(multiome_id, spatial_id, x, y), by = c('multiome_id')) %>%
    column_to_rownames(var = 'rowname')

  return(multiome.subset)

}

#' @title Convert Seurat to Anndata
#' @description Processing function to convert a Seurat object into an Anndata
#'
#' @param seurat A Seurat object. Possibly output from `MultiometoSpatialObject()`
#' @param assays An assay to convert to Anndata
#' @param dimreducs Only keep a subset of DimReducs specified here (if NULL, remove all DimReducs)
#' @param filename name of the file to save Anndata
#'
#' @return A h5ad Anndata
#' @export
ConvertAnndata <- function(seurat, assays = c('RNA'),  dimreducs = c("umap"),  filename = 'spamscOutput'){

  VariableFeatures(seurat) <- NULL
  seurat@active.assay <- assays
  if(!is.null(seurat@assays[[assays]]$scale.data)){  seurat@assays[[assays]]$scale.data <- NULL}
  seurat <- DietSeurat(seurat, assays = assays, counts = T, data = T, scale.data = F, dimreducs = dimreducs)
  SaveH5Seurat(seurat, filename = paste0(filename, '.h5Seurat'))
  Convert(paste0(filename, '.h5Seurat'), dest = "h5ad")

}



