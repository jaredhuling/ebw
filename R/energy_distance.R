

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
energy_balance <- function(trt,
                           x,
                           solver = "cccp",
                           type = c("ATE", "ATE.3", "ATT"),
                           constr.sum = FALSE,
                           quiet = TRUE,
                           alpha = NULL)
{

    type <- match.arg(type)

    trt <- as.vector(trt)

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

        res <- solvecop(lcqp, solver = solver, quiet = quiet)

        wts_quadopt <- res$x

        wts_quadopt1 <- wts_quadopt[trt == 1]
        wts_quadopt0 <- wts_quadopt[trt != 1]
        x1 <- x[trt == 1,]
        x0 <- x[trt == 0,]
        wtsall <- c(wts_quadopt0, wts_quadopt1)

        ### the overall objective function value
        final_value <- energy.dist.twoway(wtsall, x0, x1, gamma = -1)

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

        res <- solvecop(lcqp, solver = solver, quiet = quiet)

        wts_quadopt <- res$x

        wts_quadopt1 <- wts_quadopt[trt == 1]
        wts_quadopt0 <- wts_quadopt[trt != 1]
        x1 <- x[trt == 1,]
        x0 <- x[trt == 0,]
        wtsall <- c(wts_quadopt0, wts_quadopt1)

        ### the overall objective function value
        final_value <- energy.dist.threeway(wtsall, x0, x1, gamma = -1)
    }


    wts0 <- unname(wts_quadopt)[trt != 1]
    wts1 <- unname(wts_quadopt)[trt == 1]

    list(wts = unname(wts_quadopt),
         wts0 = wts0,
         wts1 = wts1,
         opt = res,
         value = final_value)
}


energy_balance2 <- function(trt,
                           x,
                           solver = "cccp",
                           type = c("ATE", "ATE.3", "ATT"),
                           constr.sum = FALSE,
                           quiet = TRUE,
                           alpha = NULL)
{

    type <- match.arg(type)

    trt <- as.factor(as.vector(trt))

    trt.levels <- levels(trt)
    K     <- length(trt.levels)
    n.vec <- unname(table(trt))
    nn    <- sum(n.vec)

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

    if (type == "two_way")
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

        res <- solvecop(lcqp, solver = solver, quiet = quiet)

        wts_quadopt <- res$x

        wts_quadopt1 <- wts_quadopt[trt == 1]
        wts_quadopt0 <- wts_quadopt[trt != 1]
        x1 <- x[trt == 1,]
        x0 <- x[trt == 0,]
        wtsall <- c(wts_quadopt0, wts_quadopt1)

        ### the overall objective function value
        final_value <- energy.dist.twoway(wtsall, x0, x1, gamma = -1)

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

        res <- solvecop(lcqp, solver = solver, quiet = quiet)

        wts_quadopt <- res$x

        wts_quadopt1 <- wts_quadopt[trt == 1]
        wts_quadopt0 <- wts_quadopt[trt != 1]
        x1 <- x[trt == 1,]
        x0 <- x[trt == 0,]
        wtsall <- c(wts_quadopt0, wts_quadopt1)

        ### the overall objective function value
        final_value <- energy.dist.threeway(wtsall, x0, x1, gamma = -1)
    }


    wts0 <- unname(wts_quadopt)[trt != 1]
    wts1 <- unname(wts_quadopt)[trt == 1]

    list(wts = unname(wts_quadopt),
         wts0 = wts0,
         wts1 = wts1,
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
        x.list[[k]]   <- x[trt_idx,]
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

        e.dist.vec[k] <- ( -2 * rbf.dist(x.list[[k]], x.target, wts.list[k], rep(1, n), gamma = gamma) +
                               rbf.dist(x.list[[k]], x.list[[k]], wts.list[k], wts.list[k], gamma = gamma) +
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
        x.list[[k]]   <- x[trt_idx,]
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

        e.dist.vec[k] <- ( -2 * rbf.dist(x.list[[k]], x.target, wts.list[k], rep(1, n), gamma = gamma) +
                               rbf.dist(x.list[[k]], x.list[[k]], wts.list[k], wts.list[k], gamma = gamma) +
                               rbf.dist(x.target, x.target, rep(1, n), rep(1, n), gamma = gamma)
        ) / gamma[1]
    }


    e.dist.pairs <- matrix(0, ncol = K, nrow = K)
    npairs <- sum(upper.tri(e.dist.pairs))

    for (k in 1:(K-1))
    {
        for (j in (k+1):K)
        {
            e.dist.pairs[k,j] <- ( -2 * rbf.dist(x.list[[k]], x.list[[j]], wts.list[k], wts.list[j], gamma = gamma) +
                                        rbf.dist(x.list[[k]], x.list[[k]], wts.list[k], wts.list[k], gamma = gamma) +
                                        rbf.dist(x.list[[j]], x.list[[j]], wts.list[j], wts.list[j], gamma = gamma)
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

