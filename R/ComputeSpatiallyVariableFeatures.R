#' @title Compute RMetric
#' @description Computes the metric at a given r (radius) value and stores
#' in meta.features from Seurat preprocessing functions
#'
#' @param mv Results of running markvario
#' @param r.metric r value at which to report the "trans" value of the markvariogram
#'
#' @return Returns a data.frame with r.metric values
#'
#' @export
ComputeRMetric <- function(mv, r.metric = 5) {
  r.metric.results <- unlist(x = lapply(
    X = mv,
    FUN = function(x) {
      x$trans[which.min(x = abs(x = x$r - r.metric))]
    }
  ))
  r.metric.results <- as.data.frame(x = r.metric.results)
  colnames(r.metric.results) <- paste0("Markvariogram_r.", r.metric)
  return(r.metric.results)
}

#' @title Compute SpatiallyV ariable Features
#' @description Compute spatially variable features from assay by
#' mark variogram and Moran's I computation
#'
#' @param data gene x cell matrix containing expression values
#' @param coords cell x 2 matrix containing spatial coordinate of cells
#' @param r.metric r value at which to report the "trans" value of the markvariogram
#'
#' @return A dataframe showing mark variogram values, Moran's I values, and
#' gene ranks based on their spatial variability.
#' @export
ComputeSpatiallyVariableFeatures <- function(data, coords, r.metric = 5 ){
  markvariogram <- FindSpatiallyVariableFeatures(
    object = data,
    spatial.location = coords,
    selection.method = c("markvariogram"),
    r.metric = r.metric
  )

  moransi <- FindSpatiallyVariableFeatures(
    object = data,
    spatial.location = coords,
    selection.method = c("moransi"),
    r.metric = r.metric
  )

  markvariogram.r <- ComputeRMetric(markvariogram)
  markvariogram.r <- markvariogram.r[order(markvariogram.r[, 1]), , drop = FALSE]
  markvariogram.r[['Markvariogram_rank']] <- 1:nrow(markvariogram.r)

  colnames(x = moransi) <- paste0("MoransI_", colnames(x = moransi))
  moransi <- moransi[order(moransi[, 2], -abs(moransi[, 1])), , drop = FALSE]
  moransi[['Moransi_rank']] <- 1:nrow(moransi)

  outs <- markvariogram.r %>% rownames_to_column(var = c('Feature')) %>%
    inner_join(moransi %>% rownames_to_column(var = c('Feature')), by = 'Feature') %>%
    mutate(Rank = rank(Markvariogram_rank + Moransi_rank)) %>%
    arrange(Rank)

  return(outs)

}


#' @title Plot Spatial Variable Ranks
#' @description Plot rank of each gene by their mark variogram and moran's I ranks.
#'
#' @param outs output dataframe from `compute_spatially_variable_features()`
#' @param rank.threshold A rank threshold to plot genes. Default:100
#' i.e. showing genes having ranks equal or less than rank.threshold.
#' rank.threshold = NA plots all genes
#' @param text.max.overlap Exclude text labels when they overlap too many other things.
#' For each text label, we count how many other text labels or other data points it overlaps,
#' and exclude the text label if it has too many overlaps. Defaults to 10.
#' @return A scatterplot showing mark variogram rank on x-axis and
#' moran's I rank on y-axis
#' @export
PlotSpatialVariableRanks <- function(outs, rank.threshold = 100, text.max.overlap = 10){

  if(is.na(rank.threshold )){
    p <- outs %>%
      ggplot(aes(x=Markvariogram_rank, y=Moransi_rank, label = Feature)) +
      geom_point() + geom_text_repel(hjust=0, vjust=0, max.overlaps = text.max.overlap, segment.color = 'grey') +
      #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
      theme_classic() + scale_x_continuous(trans='log10') +
      xlab('Markvariogram rank') + ylab('Moransi_rank')
  }else{
    p <- outs %>% filter(Rank <= rank.threshold) %>%
      ggplot(aes(x=Markvariogram_rank, y=Moransi_rank, label = Feature)) +
      geom_point() + geom_text_repel(hjust=0, vjust=0, max.overlaps = text.max.overlap, segment.color = 'grey') +
      #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
      theme_classic() + scale_x_continuous(trans='log10') +
      xlab('Markvariogram rank') + ylab('Moransi_rank')
  }

  print(p)

}

#' @title Plot Moransi Across Modality
#' @description Plot moran's I values from two assays
#'
#' @param outs output dataframe from `compute_spatially_variable_features()`
#' @param plot.moransI A moran's I value threshold to plot genes e.g. 0.3.
#' Higher the moran's I value, more spatially variable. Default: NA
#' @param top.features A number of top spatially variable genes to plot. Default:100
#' If both plot.moransI = NA & top.features = NA, plot all genes.
#' @param text.max.overlap Exclude text labels when they overlap too many other things.
#' For each text label, we count how many other text labels or other data points it overlaps,
#' and exclude the text label if it has too many overlaps. Defaults to 10.
#' @return A scatterplot showing moran's I value from assay 1 in x-axis and
#' assay 2 in y-axis
#' @export
PlotMoransiAcrossModality <- function(outs, plot.moransI = NA, top.features = 100, text.max.overlap = 10){

  # plot by moranSI
  outs.m <- outs %>% dplyr::select(Feature, MoransI_observed) %>% separate(Feature, sep = '_', into = c('Modal', 'Feature'))
  outs.m <- split(outs.m, outs.m$Modal)
  colnames(outs.m[[1]])[3] <- 'Modal1.MoransI'
  colnames(outs.m[[2]])[3] <- 'Modal2.MoransI'

  rank1 <- paste0(names(outs.m)[1], ' MoransI')
  rank2 <- paste0(names(outs.m)[2], ' MoransI')

  ranks.m <- outs.m[[1]][, c(2:3)] %>% inner_join(outs.m[[2]][, c(2:3)], by = c('Feature'))

  if(!is.na(plot.moransI)){
    p <- ranks.m %>% filter(Modal1.MoransI >= plot.moransI | Modal2.MoransI >= plot.moransI) %>%
      ggplot(aes(x=Modal1.MoransI, y=Modal2.MoransI, label = Feature)) +
      geom_point() + geom_text_repel(hjust=0, vjust=0, max.overlaps = text.max.overlap, segment.color = 'grey') +
      #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
      stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+ theme_classic() +
      xlab(rank1) + ylab(rank2)

    print(p)
  }


  if(!is.na(top.features)){
    features1 <- ranks.m %>% arrange(-Modal1.MoransI) %>% head(top.features) %>% pull(Feature)
    features2 <- ranks.m %>% arrange(-Modal2.MoransI) %>% head(top.features) %>% pull(Feature)
    p <- ranks.m %>% filter(Feature %in% unique(c(features1, features2))) %>%
      ggplot(aes(x=Modal1.MoransI, y=Modal2.MoransI, label = Feature)) +
      geom_point() + geom_text_repel(hjust=0, vjust=0, max.overlaps = text.max.overlap, segment.color = 'grey') +
      #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
      stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
      theme_classic() +
      xlab(rank1) + ylab(rank2)

    print(p)
  }


  if(is.na(plot.moransI) & is.na(top.features)){
    p <- ranks.m %>%
      ggplot(aes(x=Modal1.MoransI, y=Modal2.MoransI, label = Feature)) +
      geom_point() + geom_text_repel(hjust=0, vjust=0) +
      #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
      stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
      theme_classic() + scale_x_continuous(trans='log10') +
      xlab(rank1) + ylab(rank2)

    print(p)
  }


}

#' @title Plot Spatial Rank Across Modality
#' @description Plot spatially variable rank from two assays
#'
#' @param outs output dataframe from `compute_spatially_variable_features()`
#' @param plot.rank A rank threshold to plot genes e.g. 100.
#' Lower the rank, more spatially variable. Default: NA
#' @param top.features A number of top spatially variable genes to plot. Default:100
#' If both plot.rank = NA & top.features = NA, plot all genes.
#'
#' @return A scatterplot showing ranks from assay 1 in x-axis and
#' assay 2 in y-axis
#' @export
PlotSpatialRankAcrossModality <- function(outs, plot.rank = NA, top.features = 100){

  # plot by moranSI
  outs.m <- outs %>% dplyr::select(Feature, Rank) %>% separate(Feature, sep = '_', into = c('Modal', 'Feature'))
  outs.m <- split(outs.m, outs.m$Modal)
  colnames(outs.m[[1]])[3] <- 'Modal1.rank'
  colnames(outs.m[[2]])[3] <- 'Modal2.rank'

  rank1 <- paste0(names(outs.m)[1], ' Rank')
  rank2 <- paste0(names(outs.m)[2], ' Rank')

  ranks.m <- outs.m[[1]][, c(2:3)] %>% inner_join(outs.m[[2]][, c(2:3)], by = c('Feature'))
  ranks.m <- ranks.m %>% mutate(Modal1.rank = rank(Modal1.rank),
                                Modal2.rank = rank(Modal2.rank))

  if(!is.na(plot.rank)){
    p <- ranks.m %>% filter(Modal1.rank <= plot.rank | Modal2.rank <= plot.rank) %>%
      ggplot(aes(x=Modal1.rank, y=Modal2.rank, label = Feature)) +
      geom_point() + geom_text_repel(hjust=0, vjust=0) +
      #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
      stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
      theme_classic() + scale_x_continuous(trans='log10') +
      xlab(rank1) + ylab(rank2)


    print(p)
  }


  if(!is.na(top.features)){
    features1 <- ranks.m %>% arrange(-Modal1.rank) %>% head(top.features) %>% pull(Feature)
    features2 <- ranks.m %>% arrange(-Modal2.rank) %>% head(top.features) %>% pull(Feature)
    p <- ranks.m %>% filter(Feature %in% unique(c(features1, features2))) %>%
      ggplot(aes(x=Modal1.rank, y=Modal2.rank, label = Feature)) +
      geom_point() + geom_text_repel(hjust=0, vjust=0) +
      #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
      stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
      theme_classic() + scale_x_continuous(trans='log10') +
      xlab(rank1) + ylab(rank2)

    print(p)
  }


  if(is.na(plot.moransI) & is.na(top.features)){
    p <- ranks.m %>%
      ggplot(aes(x=Modal1.rank, y=Modal2.rank, label = Feature)) +
      geom_point() + geom_text_repel(hjust=0, vjust=0) +
      #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
      stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
      theme_classic() + scale_x_continuous(trans='log10') +
      xlab(rank1) + ylab(rank2)

    print(p)
  }


}

