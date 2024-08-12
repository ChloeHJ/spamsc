#' spamsc: A package for a comprehensive workflow designed to project multi-modal
#' single-cell data onto high-resolution spatial transcriptomic data and identify spatial features.
#'
#' The spamsc package contains functions to
#'
#' It addresses the following needs:
#' \itemize{ \item Check if multi-modal single-cell and spatial data can be integrated;
#' \item Spatially project multi-modal data;
#' \item Diagnose and evaluate projection;
#' \item Identify spatial features;
#' }
#'
#' @import Seurat
#' @import SeuratObject
#' @rawNamespace import(Matrix, except = c(head, tail))
#' @import dbscan
#' @import tibble
#' @import dplyr
#' @import ggplot2
#' @import ggpubr
#' @import viridis
#' @import reshape2
#' @import ggrepel
#' @import sf
#' @import spdep
#' @import sfdep
#' @import utils
#' @importFrom grDevices cairo_pdf dev.off
#' @importFrom stats cor fisher.test lm p.adjust
#' @importFrom magrittr %>%
#' @importFrom pheatmap pheatmap
#' @importFrom tidyr separate
#' @importFrom harmony RunHarmony
#' @rawNamespace import(matrixStats, except = count)
#'
#' @name spamsc
NULL
