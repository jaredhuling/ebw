
#' Construction of energy balancing weights
#'
#' @description Constructs energy balancing weights for estimation of average treatment effects
#'
#' @param trt vector indicating treatment assignment. Can be a factor or a vector of integers, with each
#' unique value indicating a different treatment option. When \code{method = 'ATT'}, \code{trt} MUST take
#' the value \code{1} to indicate treatment.
#' @param x matrix of covariates with number of rows equal to the length of \code{trt} and each column is a
#' \strong{pre-treatment} covariate to be balanced between treatment groups.
#' @param method What method/estimand should be used? \code{"ATE.3"} estimates weights for estimation of the
#' ATE where the energy distance criterion to be optimized involves all pairwise energy distances. This
#' method is generally better than \code{"ATE"}, which also estimates the ATE, but only involves energy distances
#' bbetween each treatment group and the full sample. \code{"ATT"} estimates weights for estimating the ATT
#' @param standardized logical scalar which results in standardization of the covariates when it takes the value \code{TRUE}
#' and when it takes the value \code{FALSE}, covariates are not standardized. The default is \code{TRUE}, which is
#' highly recommended
#' @param verbose should we print out intermediate results of optimization process? \code{TRUE} or \code{FALSE}
#' @param constr.sum should each weight be constrained to be less than \code{10 * nrow(x) ^ (1/3)}? Defaults to \code{FALSE}.
#' This option is generally not needed.
#' @return An object of class \code{"energy_balancing_weights"} with elements:
#' \item{weights}{A vector of length \code{nrow(x)} containing the estimates sample weights }
#' \item{trt}{Treatment vector}
#' \item{estimand}{The estimand requested}
#' \item{method}{The method used}
#' \item{covs}{The covariates that were used for balancing}
#' \item{energy_dist_unweighted}{The energy distance of the raw data (ie with all weights = 1)}
#' \item{energy_dist_optimized}{The weighted energy distance using the optimal energy balancing weights}
#' \item{opt}{The optimization object returned by \code{solvecop()}}
#'
#' @export
energy_balance <- function(trt,
                           x,
                           method      = c("ATE.3", "ATE", "ATT", "overlap"),
                           standardize = TRUE,
                           verbose     = FALSE,
                           constr.sum  = FALSE)
{

    type <- match.arg(method)

    if (standardize)
    {
        x.orig <- x
        x <- scale(x)
    }

    trt <- as.factor(as.vector(trt))

    trt.levels <- levels(trt)
    K          <- length(trt.levels)
    n.vec      <- unname(table(trt))
    nn         <- sum(n.vec)

    solver = "cccp"

    if (K > 2 & type == "ATT")
    {
        stop("type == 'ATT' cannot be used for treatments with multiple levels")
    }

    if (K > 2 & type == "overlap")
    {
        stop("type == 'overlap' cannot be used for treatments with multiple levels")
    }

    if (type == "ATT")
    {
        if (trt.levels[2] != "1")
        {
            stop("for type == 'ATT', 'trt' must take the value 1 for treated units ")
        }
    }

    alpha       = NULL
    QQ_all <- rdist(x)

    if (!is.null(alpha))
    {
        QQ_all <- -exp(alpha * QQ_all ^ 2)
    }

    N  <- nn + sum(n.vec)


    AA.constr.mat <- matrix(0, nrow = 1, ncol = nn)

    AA.list <- rep(list(AA.constr.mat), K)

    for (k in 1:K)
    {
        AA.list[[k]][trt == trt.levels[k]] <- 1
    }

    if (type == "ATE")
    {

        ## two_way balances the treatment arm to the full population
        ## and              the control   arm to the full population
        for (k in 1:K)
        {
            if (k == 1)
            {
                trt_ind <- 1 * (trt == trt.levels[k])
                QQ <- -trt_ind * t( trt_ind * t(QQ_all)) / n.vec[k] ^ 2

                aa <- 2 * as.vector(rowSums(trt_ind * QQ_all)) / (n.vec[k] * nn)
            } else
            {
                trt_ind <- 1 * (trt == trt.levels[k])
                QQ <- QQ - trt_ind * t( trt_ind * t(QQ_all)) / n.vec[k] ^ 2

                aa <- aa + 2 * as.vector(rowSums(trt_ind * QQ_all)) / (n.vec[k] * nn)
            }
        }



        rownames(QQ) <- paste(1:NROW(QQ))

        qf <- quadfun(Q = QQ, a = aa, id = rownames(QQ)) #quadratic obj.
        lb <- lbcon(0, id = rownames(QQ)) #lower bound

        AA           <- do.call(rbind, AA.list)
        rownames(AA) <- paste0("eq", 1:nrow(AA))
        sum.constr   <- n.vec


        lc <- lincon(A   = AA,
                     dir = rep("==", length(sum.constr)),
                     val = sum.constr,
                     id  = rownames(QQ)) #linear constraint


        if (constr.sum)
        {
            ub <- ubcon(10 * nn ^ (1/3), id = rownames(QQ))
            lcqp <- cop( f = qf, lb = lb, lc = lc, ub = ub)
        } else
        {
            lcqp <- cop( f = qf, lb = lb, lc = lc)
        }

        res <- solvecop(lcqp, solver = solver, quiet = !verbose)

        wts_quadopt <- unname(res$x)

        ### the overall objective function value
        final_value <- energy.dist.1bal.trt(unname(res$x), x, trt, gamma = -1)

        ### the unweighted objective function value
        unweighted_value <- energy.dist.1bal.trt(rep(1, length(res$x)), x, trt, gamma = -1)

    } else if (type == "ATE.3")
    {
        ## three_way balances the treatment arm to the full population
        ##                    the control   arm to the full population
        ## and                the treatment arm to the control arm


        for (k in 1:K)
        {
            if (k == 1)
            {
                trt_ind <- 1 * (trt == trt.levels[k])
                QQ <- -trt_ind * t( trt_ind * t(QQ_all)) / n.vec[k] ^ 2

                aa <- 2 * as.vector(rowSums(trt_ind * QQ_all)) / (n.vec[k] * nn)
            } else
            {
                trt_ind <- 1 * (trt == trt.levels[k])
                QQ <- QQ - trt_ind * t( trt_ind * t(QQ_all)) / n.vec[k] ^ 2

                aa <- aa + 2 * as.vector(rowSums(trt_ind * QQ_all)) / (n.vec[k] * nn)
            }
        }

        for (k in 1:(K - 1))
        {
            trt_ind_k <- 1 * (trt == trt.levels[k])
            for (j in (k + 1):K)
            {
                trt_ind_j <- 1 * (trt == trt.levels[j])
                QQ <- QQ + 2 * (trt_ind_k) * t( trt_ind_j * t(QQ_all)) / (n.vec[k] * n.vec[j])

                QQ <- QQ - trt_ind_j * t( trt_ind_j * t(QQ_all)) / n.vec[j] ^ 2
                QQ <- QQ - trt_ind_k * t( trt_ind_k * t(QQ_all)) / n.vec[k] ^ 2
            }
        }


        AA           <- do.call(rbind, AA.list)
        rownames(AA) <- paste0("eq", 1:nrow(AA))
        sum.constr   <- n.vec

        rownames(QQ) <- paste(1:NROW(QQ))

        qf <- quadfun(Q = QQ, a = aa, id = rownames(QQ)) #quadratic obj.
        lb <- lbcon(0, id = rownames(QQ)) #lower bound


        lc <- lincon(A   = AA,
                     dir = rep("==", length(sum.constr)),
                     val = sum.constr,
                     id  = rownames(QQ)) #linear constraint


        if (constr.sum)
        {
            ub <- ubcon(10 * nn ^ (1/3), id = rownames(QQ))
            lcqp <- cop( f = qf, lb = lb, lc = lc, ub = ub)
        } else
        {
            lcqp <- cop( f = qf, lb = lb, lc = lc)
        }

        res <- solvecop(lcqp, solver = solver, quiet = !verbose)

        wts_quadopt <- unname(res$x)

        ### the overall objective function value
        final_value <- energy.dist.pairbal.trt(unname(res$x), x, trt, gamma = -1)

        ### the unweighted objective function value
        unweighted_value <- energy.dist.pairbal.trt(rep(1, length(res$x)), x, trt, gamma = -1)
    } else if (type == "ATT")
    {
        ## these are the untreated patients
        trt_ind <- 1 * (trt == trt.levels[1])
        QQ <- -trt_ind * t( trt_ind * t(QQ_all)) / n.vec[1] ^ 2

        aa <- 2 * as.vector(rowSums(trt_ind * QQ_all)) / (n.vec[1] * nn)
    } else if (type == "overlap")
    {
        ## these are the treated patients
        trt_ind <- 1 * (trt == trt.levels[2])

        QQ1 <- trt_ind * t( trt_ind * t(QQ_all)) / n.vec[2] ^ 2
        QQ0 <- (1 - trt_ind) * t( (1 - trt_ind) * t(QQ_all)) / n.vec[1] ^ 2



        QQ_both <- (trt_ind) * t( (1 - trt_ind) * t(QQ_all)) / (n.vec[1] * n.vec[2])

        QQ <- 2 * QQ_both - (QQ1 + QQ0)

        AA <- matrix(1, nrow = 1, ncol = nn)
        rownames(AA) <- "eq"
        rownames(QQ) <- paste(1:NROW(QQ))


        qf <- quadfun(Q = QQ, id = rownames(QQ)) #quadratic obj.
        lb <- lbcon(0, id = rownames(QQ))
        ub <- ubcon(1, id = rownames(QQ))


        lcqp <- cop( f = qf, lb = lb, ub = ub)

        res <- solvecop(lcqp, solver = solver, quiet = !verbose)

        wts_quadopt <- unname(res$x)

        #wts_quadopt[trt_ind == 1] <- 1


        ### the overall objective function value
        final_value <- weighted.energy.dist.trt2ctrl(unname(res$x), x, trt, gamma = -1)

        ### the unweighted objective function value
        unweighted_value <- weighted.energy.dist.trt2ctrl(rep(1, length(res$x)), x, trt, gamma = -1)
    }


    if (type %in% c("ATE", "ATE.3"))
    {
        estimand <- "ATE"
    } else if (type %in% c("ATT"))
    {
        estimand <- "ATT"
    } else if (type == "overlap")
    {
        estimand <- "overlap"
    }


    ret <- list(weights                = wts_quadopt,
                treat                  = trt,
                estimand               = estimand,
                method                 = type,
                covs                   = if (standardize) {x.orig} else {x},
                energy_dist_unweighted = unweighted_value,
                energy_dist_optimized  = final_value,
                opt                    = res
                )
    class(ret) <- c("energy_balancing_weights")
}



# energy.dist.1bal.trt <- function(wts, x, trt, gamma = -1, normalize.wts = FALSE)


#' @export
weighted_energy_distance <- function(weights = rep(1, NROW(x)),
                                     x,
                                     trt,
                                     type = c("two_way", "three_way"),
                                     normalize.wts = TRUE)
{
    type <- match.arg(type)

    stopifnot(length(weights) == NROW(x))
    stopifnot(length(weights) == NROW(trt))

    if (type == "two_way")
    {
        return(energy.dist.1bal.trt(wts = weights, x = x, trt = trt, gamma = -1,
                                    normalize.wts = normalize.wts))
    } else
    {
        return(energy.dist.pairbal.trt(wts = weights, x = x, trt = trt, gamma = -1,
                                       normalize.wts = normalize.wts))
    }
}






##  This is the main function for estimating
##  "energy balancing weights"
##
##
## trt is a treatment vector taking values 0 and 1
## x is a matrix of covariates (number of rows = length of trt)
## solver is the solver type passed to solvecop()
## type: should we only balance each arm to the population? (two_way)
##       or balance that plus trt to control? (three_way)
## quiet: should we tell solvecop() to be quiet?
energy_balance_old <- function(trt,
                               x,
                               solver = "cccp",
                               type = c("ATE", "ATE.3", "ATT"),
                               constr.sum = FALSE,
                               standardize = TRUE,
                               verbose = FALSE,
                               alpha = NULL)
{

    type <- match.arg(type)

    trt <- as.vector(trt)


    if (standardize)
    {
        x <- scale(x)
    }

    QQ_all <- rdist(x)

    if (!is.null(alpha))
    {
        QQ_all <- -exp(alpha * QQ_all ^ 2)
    }

    n1 <- sum(trt)
    n0 <- sum(1-trt)
    nn <- length(trt)

    N  <- nn + n0 + n1


    AA1 <- matrix(1, nrow = 1, ncol = nn)
    AA0 <- matrix(1, nrow = 1, ncol = nn)

    AA1[trt != 1] <- 0
    AA0[trt == 1] <- 0

    if (type == "two_way")
    {
        ## two_way balances the treatment arm to the full population
        ## and              the control   arm to the full population
        QQ1 <- trt * t( trt * t(QQ_all)) / n1 ^ 2
        QQ0 <- (1 - trt) * t( (1 - trt) * t(QQ_all)) / n0 ^ 2

        aa1 <- 2 * as.vector(rowSums(trt * QQ_all)) / (n1 * nn)
        aa0 <- 2 * as.vector(rowSums((1 - trt) * QQ_all)) / (n0 * nn)


        aa <- aa1 + aa0
        QQ <- -(QQ1 + QQ0)
        rownames(QQ) <- paste(1:NROW(QQ))

        qf <- quadfun(Q = QQ, a = aa, id = rownames(QQ)) #quadratic obj.
        lb <- lbcon(0, id = rownames(QQ)) #lower bound

        # if (sum.to.one == "all")
        # {
        #     AA <- matrix(1, nrow = 1, ncol = nn)
        #     rownames(AA) <- "eq"
        #     sum.constr <- nn
        # } else
        # {
        AA <- rbind(AA1, AA0)
        rownames(AA) <- paste0("eq", 1:nrow(AA))
        sum.constr <- c(n1, n0)
        # }


        lc <- lincon(A   = AA,
                     dir = rep("==", length(sum.constr)),
                     val = sum.constr,
                     id  = rownames(QQ)) #linear constraint



        if (constr.sum)
        {
            ub <- ubcon(10 * nn ^ (1/3), id = rownames(QQ))
            lcqp <- cop( f = qf, lb = lb, lc = lc, ub = ub)
        } else
        {
            lcqp <- cop( f = qf, lb = lb, lc = lc)
        }

        res <- solvecop(lcqp, solver = solver, quiet = !verbose)

        wts_quadopt <- res$x

        ### the overall objective function value
        final_value <- energy.dist.1bal.trt(res$x, x, trt, gamma = -1)

        ### the unweighted objective function value
        unweighted_value <- energy.dist.1bal.trt(rep(1, length(res$x)), x, trt, gamma = -1)

    } else
    {
        ## three_way balances the treatment arm to the full population
        ##                    the control   arm to the full population
        ## and                the treatment arm to the control arm
        QQ1 <- trt * t( trt * t(QQ_all)) / n1 ^ 2
        QQ0 <- (1 - trt) * t( (1 - trt) * t(QQ_all)) / n0 ^ 2


        QQ_both <- (trt) * t( (1 - trt) * t(QQ_all)) / (n1 * n0)

        aa1 <- 2 * as.vector(rowSums(trt * QQ_all)) / (n1 * nn)
        aa0 <- 2 * as.vector(rowSums((1 - trt) * QQ_all)) / (n0 * nn)


        aa <- aa1 + aa0
        QQ <- 2 * QQ_both - 2 * (QQ1 + QQ0)

        AA <- matrix(1, nrow = 1, ncol = nn)
        rownames(AA) <- "eq"
        rownames(QQ) <- paste(1:NROW(QQ))


        qf <- quadfun(Q = QQ, a = aa, id = rownames(QQ)) #quadratic obj.
        lb <- lbcon(0, id = rownames(QQ)) #lower bound

        # if (sum.to.one == "all")
        # {
        #     AA <- matrix(1, nrow = 1, ncol = nn)
        #     rownames(AA) <- "eq"
        #     sum.constr <- nn
        # } else
        # {
        AA <- rbind(AA1, AA0)
        rownames(AA) <- paste0("eq", 1:nrow(AA))
        sum.constr <- c(n1, n0)
        # }


        lc <- lincon(A   = AA,
                     dir = rep("==", length(sum.constr)),
                     val = sum.constr,
                     id  = rownames(QQ)) #linear constraint


        if (constr.sum)
        {
            ub <- ubcon(10 * nn ^ (1/3), id = rownames(QQ))
            lcqp <- cop( f = qf, lb = lb, lc = lc, ub = ub)
        } else
        {
            lcqp <- cop( f = qf, lb = lb, lc = lc)
        }

        res <- solvecop(lcqp, solver = solver, quiet = !verbose)

        wts_quadopt <- res$x

        ### the overall objective function value
        final_value <- energy.dist.pairbal.trt(res$x, x, trt, gamma = -1)

        ### the unweighted objective function value
        unweighted_value <- energy.dist.pairbal.trt(rep(1, length(res$x)), x, trt, gamma = -1)
    }

    list(wts = unname(res$x),
         energy_dist_unweighted = unweighted_value,
         energy_dist_optimized  = final_value,
         opt = res,
         value = final_value)
}



## calculates weighted distance
## when gamma is negative, it calculates euclidean distance,
## otherwise it uses a gaussian kernel
rbf.dist <- function(x, y, weights.x = NULL, weights.y = NULL, gamma = -1)
{
    x   <- as.matrix(x)
    y   <- as.matrix(y)
    n.x <- NROW(x)
    n.y <- NROW(y)

    if (is.null(weights.x))
    {
        weights.x <- rep(1, n.x)
    }
    if (is.null(weights.y))
    {
        weights.y <- rep(1, n.y)
    }

    weights.x <- weights.x / mean(weights.x)
    weights.y <- weights.y / mean(weights.y)

    if (gamma < 0)
    {
        rd <- rdist(x, y)
    } else
    {
        rd <- exp(-gamma * rdist(x, y) ^ 2)
    }

    sum(weights.x * t(weights.y * t(rd))) / (n.x * n.y)
}


energy.dist <- function(x0, x1, gamma = -1)
{
    n.x0 <- NROW(x0)
    n.x1 <- NROW(x1)

    wts.x0 <- NULL
    wts.x1 <- NULL

    ( -2 *  rbf.dist(x0, x1, wts.x0, wts.x1, gamma = gamma) +
            rbf.dist(x0, x0, wts.x0, wts.x0, gamma = gamma) +
            rbf.dist(x1, x1, wts.x1, wts.x1, gamma = gamma)
    ) / gamma[1]
}


weighted.energy.dist <- function(wts, x0, x1, gamma = 1)
{
    n.x0 <- NROW(x0)
    n.x1 <- NROW(x1)

    wts.x0 <- wts[1:n.x0]
    wts.x1 <- wts[-(1:n.x0)]

    ( -2 *  rbf.dist(x0, x1, wts.x0, wts.x1, gamma = gamma) +
            rbf.dist(x0, x0, wts.x0, wts.x0, gamma = gamma) +
            rbf.dist(x1, x1, wts.x1, wts.x1, gamma = gamma)
    ) / gamma[1]
}

weighted.energy.dist.trt2ctrl <- function(wts, x, trt, gamma = 1, normalize.wts = FALSE)
{

    trt_idx <- trt == 1

    x1 <- x[trt_idx,,drop=FALSE]
    x0 <- x[!trt_idx,,drop=FALSE]

    n.0  <- NROW(x0)
    n.1  <- NROW(x1)

    if (normalize.wts)
    {
        wts.x0 <- wts[!trt_idx] / mean(wts[!trt_idx])
        wts.x1 <- wts[trt_idx] / mean(wts[trt_idx])
    } else
    {
        wts.x0 <- wts[!trt_idx]
        wts.x1 <- wts[trt_idx]
    }

    ( -2 *  rbf.dist(x0, x1, wts.x0, wts.x1, gamma = gamma) +
            rbf.dist(x0, x0, wts.x0, wts.x0, gamma = gamma) +
            rbf.dist(x1, x1, wts.x1, wts.x1, gamma = gamma)
    ) / gamma[1]
}


weighted.energy.dist.trt2full <- function(wts, x, trt, gamma = 1, normalize.wts = FALSE)
{

    trt_idx <- trt == 1

    x1 <- x[trt_idx,,drop=FALSE]

    n    <- NROW(x)
    n.1  <- NROW(x1)

    x.all <- rbind(x1, x[!trt_idx,,drop=FALSE])

    wts.x  <- rep(1, n)
    wts.x1 <- wts[trt_idx]

    if (normalize.wts)
    {
        wts.x  <- rep(1, n)
        wts.x1 <- wts[trt_idx] / mean(wts[trt_idx])
    } else
    {
        wts.x  <- rep(1, n)
        wts.x1 <- wts[trt_idx]
    }

    ( -2 *  rbf.dist(x.all, x1, wts.x, wts.x1, gamma = gamma) +
            rbf.dist(x.all, x.all, wts.x, wts.x, gamma = gamma) +
            rbf.dist(x1, x1, wts.x1, wts.x1, gamma = gamma)
    ) / gamma[1]
}



energy.dist.1bal.trt <- function(wts, x, trt, gamma = -1, normalize.wts = FALSE)
{
    if (!is.factor(trt))
    {
        trt <- as.factor(trt)
    }

    trt.levels <- levels(trt)
    K     <- length(trt.levels)
    n.vec <- unname(table(trt))
    n    <- sum(n.vec)

    wts.list <- x.list <- vector(mode = "list", length = K)

    for (k in 1:K)
    {
        trt_idx       <- trt == trt.levels[k]
        x.list[[k]]   <- x[trt_idx,,drop=FALSE]
        if (normalize.wts)
        {
            wts.list[[k]] <- wts[trt_idx] / mean(wts[trt_idx])
        } else
        {
            wts.list[[k]] <- wts[trt_idx]
        }
    }

    n    <- sum(n.vec)
    N    <- n + sum(n.vec)

    ## full population
    x.target <- do.call(rbind, x.list)

    e.dist.vec <- numeric(k)

    for (k in 1:K)
    {
        ## dist between w*x.trt and (x.trt, x1)

        e.dist.vec[k] <- ( -2 * rbf.dist(x.list[[k]], x.target, wts.list[[k]], rep(1, n), gamma = gamma) +
                               rbf.dist(x.list[[k]], x.list[[k]], wts.list[[k]], wts.list[[k]], gamma = gamma) +
                               rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
                         ) / gamma[1]
    }

    ## dist between w*x0 and (x0, x1)

    e.dist.vec.scaled <- n * e.dist.vec * n.vec / (2 * N)

    # between sample energy statistic (6.7 of http://www.ericbunch.org/static/Szekely_estats.pdf)
    #e.dist.0 * ((n.x0 * n) / (2 * N)) + e.dist.1 * ((n.x1 * n) / (2 * N)) ## + e.dist.01 * ((n.x0 * n.x1) / (2 * N))

    sum(e.dist.vec.scaled)
}


energy.dist.pairbal.trt <- function(wts, x, trt, gamma = -1, normalize.wts = FALSE)
{
    if (!is.factor(trt))
    {
        trt <- as.factor(trt)
    }

    trt.levels <- levels(trt)
    K     <- length(trt.levels)
    n.vec <- unname(table(trt))
    n    <- sum(n.vec)

    wts.list <- x.list <- vector(mode = "list", length = K)

    for (k in 1:K)
    {
        trt_idx       <- trt == trt.levels[k]
        x.list[[k]]   <- x[trt_idx,,drop=FALSE]
        if (normalize.wts)
        {
            wts.list[[k]] <- wts[trt_idx] / mean(wts[trt_idx])
        } else
        {
            wts.list[[k]] <- wts[trt_idx]
        }
    }

    n    <- sum(n.vec)
    N    <- n + sum(n.vec)

    ## full population
    x.target <- do.call(rbind, x.list)

    e.dist.vec <- numeric(k)

    for (k in 1:K)
    {
        ## dist between w*x.trt and (x.trt, x1)

        e.dist.vec[k] <- ( -2 * rbf.dist(x.list[[k]], x.target, wts.list[[k]], rep(1, n), gamma = gamma) +
                               rbf.dist(x.list[[k]], x.list[[k]], wts.list[[k]], wts.list[[k]], gamma = gamma) +
                               rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
        ) / gamma[1]
    }


    e.dist.pairs <- matrix(0, ncol = K, nrow = K)
    npairs <- sum(upper.tri(e.dist.pairs))

    for (k in 1:(K-1))
    {
        for (j in (k+1):K)
        {
            e.dist.pairs[k,j] <- ( -2 * rbf.dist(x.list[[k]], x.list[[j]], wts.list[[k]], wts.list[[j]], gamma = gamma) +
                                        rbf.dist(x.list[[k]], x.list[[k]], wts.list[[k]], wts.list[[k]], gamma = gamma) +
                                        rbf.dist(x.list[[j]], x.list[[j]], wts.list[[j]], wts.list[[j]], gamma = gamma)
            ) / gamma[1]

            ## scale the e dist
            e.dist.pairs[k,j] <- e.dist.pairs[k,j] * n.vec[k] * n.vec[j] / (2 * N)
        }
    }


    ## dist between w*x0 and (x0, x1)

    e.dist.vec.scaled <- n * e.dist.vec * n.vec / (2 * N)

    # between sample energy statistic (6.7 of http://www.ericbunch.org/static/Szekely_estats.pdf)
    #e.dist.0 * ((n.x0 * n) / (2 * N)) + e.dist.1 * ((n.x1 * n) / (2 * N)) ## + e.dist.01 * ((n.x0 * n.x1) / (2 * N))
    ## e.dist.01 * ((n.x0 * n.x1) / (2 * N))

    sum(e.dist.vec.scaled) + sum(e.dist.pairs[upper.tri(e.dist.pairs)])
}



energy.dist.twoway <- function(wts, x0, x1, gamma = -1)
{
    n.x0 <- NROW(x0)
    n.x1 <- NROW(x1)

    n    <- n.x0 + n.x1
    N    <- n + n.x0 + n.x1

    ## full population
    x.target <- rbind(x0, x1)

    wts.x0 <- wts[1:n.x0]
    wts.x1 <- wts[-(1:n.x0)]


    ## dist between w*x0 and (x0, x1)

    e.dist.0 <- ( -2 * rbf.dist(x0, x.target, wts.x0, rep(1, n), gamma = gamma) +
                      rbf.dist(x0, x0, wts.x0, wts.x0, gamma = gamma) +
                      rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
    ) / gamma[1]#ifelse(gamma[1] < 0, 1, gamma[1])

    ## dist between w*x1 and (x0, x1)

    e.dist.1 <- ( -2 * rbf.dist(x1, x.target, wts.x1, rep(1, n), gamma = gamma) +
                      rbf.dist(x1, x1, wts.x1, wts.x1, gamma = gamma) +
                      rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
    ) / gamma[1]


    # between sample energy statistic (6.7 of http://www.ericbunch.org/static/Szekely_estats.pdf)
    e.dist.0 * ((n.x0 * n) / (2 * N)) + e.dist.1 * ((n.x1 * n) / (2 * N)) ## + e.dist.01 * ((n.x0 * n.x1) / (2 * N))
}


energy.dist.twoway.trt <- function(wts, x, trt, gamma = -1, normalize.wts = FALSE)
{
    trt_idx <- trt == 1

    x1 <- x[trt_idx,]
    x0 <- x[!trt_idx,]

    wts.x0 <- wts[!trt_idx]
    wts.x1 <- wts[trt_idx]

    if (normalize.wts)
    {
        wts.x0 <- wts.x0 / mean(wts.x0)
        wts.x1 <- wts.x1 / mean(wts.x1)
    }

    n.x0 <- NROW(x0)
    n.x1 <- NROW(x1)

    n    <- n.x0 + n.x1
    N    <- n + n.x0 + n.x1

    ## full population
    x.target <- rbind(x0, x1)


    ## dist between w*x0 and (x0, x1)

    e.dist.0 <- ( -2 * rbf.dist(x0, x.target, wts.x0, rep(1, n), gamma = gamma) +
                      rbf.dist(x0, x0, wts.x0, wts.x0, gamma = gamma) +
                      rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
    ) / gamma[1]#ifelse(gamma[1] < 0, 1, gamma[1])

    ## dist between w*x1 and (x0, x1)

    e.dist.1 <- ( -2 * rbf.dist(x1, x.target, wts.x1, rep(1, n), gamma = gamma) +
                      rbf.dist(x1, x1, wts.x1, wts.x1, gamma = gamma) +
                      rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
    ) / gamma[1]


    # between sample energy statistic (6.7 of http://www.ericbunch.org/static/Szekely_estats.pdf)
    e.dist.0 * ((n.x0 * n) / (2 * N)) + e.dist.1 * ((n.x1 * n) / (2 * N)) ## + e.dist.01 * ((n.x0 * n.x1) / (2 * N))
}



energy.dist.threeway <- function(wts, x0, x1, gamma = -1)
{
    n.x0 <- NROW(x0)
    n.x1 <- NROW(x1)

    n    <- n.x0 + n.x1

    N    <- n + n.x0 + n.x1

    x.target <- rbind(x0, x1)

    wts.x0 <- wts[1:n.x0]
    wts.x1 <- wts[-(1:n.x0)]


    ## dist between w*x0 and (x0, x1)

    e.dist.0 <- ( -2 * rbf.dist(x0, x.target, wts.x0, rep(1, n), gamma = gamma) +
                      rbf.dist(x0, x0, wts.x0, wts.x0, gamma = gamma) +
                      rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
    ) / gamma[1]

    ## dist between w*x1 and (x0, x1)

    e.dist.1 <- ( -2 * rbf.dist(x1, x.target, wts.x1, rep(1, n), gamma = gamma) +
                      rbf.dist(x1, x1, wts.x1, wts.x1, gamma = gamma) +
                      rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
    ) / gamma[1]

    ## dist between w*x0 and w*x1

    e.dist.01 <- ( -2 * rbf.dist(x0, x1, wts.x0, wts.x1, gamma = gamma) +
                       rbf.dist(x0, x0, wts.x0, wts.x0, gamma = gamma) +
                       rbf.dist(x1, x1, wts.x1, wts.x1, gamma = gamma)
    ) / gamma[1]


    # between sample energy statistic (6.7 of http://www.ericbunch.org/static/Szekely_estats.pdf)
    e.dist.0 * ((n.x0 * n) / (2 * N)) + e.dist.01 * ((n.x0 * n.x1) / (2 * N)) + e.dist.1 * ((n.x1 * n) / (2 * N))
}


energy.dist.threeway.trt <- function(wts, x, trt, gamma = -1, normalize.wts = FALSE)
{

    trt_idx <- trt == 1

    x1 <- x[trt_idx,]
    x0 <- x[!trt_idx,]

    wts.x0 <- wts[!trt_idx]
    wts.x1 <- wts[trt_idx]

    if (normalize.wts)
    {
        wts.x0 <- wts.x0 / mean(wts.x0)
        wts.x1 <- wts.x1 / mean(wts.x1)
    }

    n.x0 <- NROW(x0)
    n.x1 <- NROW(x1)

    n    <- n.x0 + n.x1

    N    <- n + n.x0 + n.x1

    x.target <- rbind(x0, x1)


    ## dist between w*x0 and (x0, x1)

    e.dist.0 <- ( -2 * rbf.dist(x0, x.target, wts.x0, rep(1, n), gamma = gamma) +
                      rbf.dist(x0, x0, wts.x0, wts.x0, gamma = gamma) +
                      rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
    ) / gamma[1]

    ## dist between w*x1 and (x0, x1)

    e.dist.1 <- ( -2 * rbf.dist(x1, x.target, wts.x1, rep(1, n), gamma = gamma) +
                      rbf.dist(x1, x1, wts.x1, wts.x1, gamma = gamma) +
                      rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
    ) / gamma[1]

    ## dist between w*x0 and w*x1

    e.dist.01 <- ( -2 * rbf.dist(x0, x1, wts.x0, wts.x1, gamma = gamma) +
                       rbf.dist(x0, x0, wts.x0, wts.x0, gamma = gamma) +
                       rbf.dist(x1, x1, wts.x1, wts.x1, gamma = gamma)
    ) / gamma[1]


    # between sample energy statistic (6.7 of http://www.ericbunch.org/static/Szekely_estats.pdf)
    e.dist.0 * ((n.x0 * n) / (2 * N)) + e.dist.01 * ((n.x0 * n.x1) / (2 * N)) + e.dist.1 * ((n.x1 * n) / (2 * N))
}

