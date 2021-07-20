
#' Calculation of Hajek estimate of the average treatment effect
#'
#' @description Calculates the Hajek weighted estimate of the average treatment effect
#'
#' @param y vector of responses
#' @param trt vector indicating treatment assignment of length equal to the length of \code{y}. It should
#' take the value 1 for treated and any other value to indicate not treated
#' @param weights a vector of sample weights of length equal to the length of \code{y}
#' @return an estimate of a causal contrast
#'
#' @export
ipw_est <- function(y, trt, weights)
{
    stopifnot(length(trt) == length(y))
    stopifnot(length(trt) == length(weights))

    if (is.factor(trt))
    {
        trt.levels <- levels(trt)
    } else
    {
        trt.levels <- unique(trt)
    }

    if (!(1 %in% trt.levels))
    {
        stop("'trt' must take value 1 to indicate treated")
    }

    sum( (trt * y * weights)[trt == 1] ) / sum(weights[trt == 1]) -
        sum( ((1 - trt) * y * weights)[trt != 1] ) / sum(weights[trt != 1])
}


## bivariate (weighted) empirical CDF functions
ecdf2d     <- function(obj, x, y) sum( obj[,1] < x & obj[,2] < y)/nrow(obj)
ecdf2d.wtd <- function(obj, x, y, wts = rep(1, NROW(obj))) sum( wts * (obj[,1] < x & obj[,2] < y)  ) / sum(wts)

## univariate (weighted) ECDF
ecdf.wtd <- function(obj, x, wts = rep(1, NROW(obj)))
{
    if (is.null(dim(obj)))
    {
        obj <- matrix(obj, ncol = 1)
    }
    nc <- NCOL(obj)
    nr <- NROW(obj)
    xmat <- matrix(rep(x, nr), ncol = nc, byrow = TRUE)
    sum( wts * (rowSums(obj < xmat) == nc)  ) / sum(wts)
}

ecdf1d.wtd <- function(obj, x, wts = rep(1, NROW(obj))) sum( wts * (obj < x)  ) / sum(wts)


## evaluate the multivariate (weighted) CDFs on a grid
bivar.wcdf.grid <- function(obj, wts, g1, g2)
{
    cdf_mat <- matrix(0, nrow = length(g1), ncol = length(g2))

    for (xc1 in 1:length(g1))
    {
        for (xc2 in 1:length(g2))
        {
            cdf_mat[xc1, xc2] <- ecdf2d.wtd(obj, x = g1[xc1], y = g2[xc2], wts = wts)
        }
    }
    cdf_mat
}

## evaluate the multivariate (weighted) CDFs on a grid
multivar.wcdf.grid <- function(obj, wts, gmat)
{
    cdf_vec <- numeric(nrow(gmat))

    for (i in 1:nrow(gmat))
    {
        cdf_vec[i] <- ecdf.wtd(obj, x = gmat[i,], wts = wts)
    }
    cdf_vec
}


## evaluate the multivariate (weighted) CDFs on a grid
univar.wcdf.grid <- function(obj, wts, g1)
{
    cdf_mat <- matrix(0, nrow = length(g1), ncol = 1)

    for (xc1 in 1:length(g1))
    {
        cdf_mat[xc1, 1] <- ecdf1d.wtd(obj, x = g1[xc1], wts = wts)
    }
    cdf_mat
}

bivar.cdf.grid <- function(obj, g1, g2)
{
    cdf_mat <- matrix(0, nrow = length(g1), ncol = length(g2))

    for (xc1 in 1:length(g1))
    {
        for (xc2 in 1:length(g2))
        {
            cdf_mat[xc1, xc2] <- ecdf2d(obj, x = g1[xc1], y = g2[xc2])
        }
    }
    cdf_mat
}

