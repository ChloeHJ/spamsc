#' @title Compute autocorrelation
#' @description Compute spatial autocorrelation of genes
#'
#' @param coords cell x 2 matrix containing spatial coordinate of cells
#' @param data gene x cell matrix containing expression values
#'
#' @return A gene x gene matrix showing spatial autocorrelation represented by
#' global bivariate Moran's index
#'
#' @export
ComputeAutocorrelation <- function(data, coords){

  nb <- knn2nb(knearneigh(coords, k = 5))
  wt <- st_weights(nb, style = "W")
  # nb <- knn2nb(knearneigh(coords, k = 5))
  # listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  mtx <- matrix(nrow=nrow(data), ncol=nrow(data), dimnames = list(rownames(data), rownames(data)))



  for( i in 1:nrow(data)){
    if (i %% 30 == 0) {
      cat(sprintf("Computing %dth genes\n", i))
    }

    for( j in 1:nrow(data)){
      tryCatch({
        bv <- global_moran_bv(data[i, ], data[j, ], nb, wt, scale = TRUE)
        mtx[i, j] <- bv$t0

        # z_var1 <- scale(data[i, ])
        # z_var2 <- scale(data[j, ])
        # lag_z_var2 <- lag.listw(listw, z_var2)
        # I <- sum(z_var1 * lag_z_var2) / sum(z_var1^2)
        # mtx[i, j] <- I
        #
      }, error=function(x){
        mtx[i, j] <- NA
      })
    }
  }
  return(mtx)

}
