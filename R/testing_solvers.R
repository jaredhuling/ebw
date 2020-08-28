
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
#' @param lambda tuning parameter for the penalty on the sum of squares of the weights
#' @param verbose should we print out intermediate results of optimization process? \code{TRUE} or \code{FALSE}
#' @param max.constr should each weight be constrained to be less than \code{10 * nrow(x) ^ (1/3)}? Defaults to \code{FALSE}.
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
#' \seealso \code{\link[ebw]{print.energy_balancing_weights}} for printing of fitted energy balancing objects
#'
#' @examples
#'
#' n <- 100
#' p <- 5
#'
#' set.seed(1)
#'
#' dat <- sim_confounded_data(n.obs = n, n.vars = p, AR.cor = 0.75,
#'                            propensity.model = "IV", y.model = "A")
#'
#' x   <- dat$x
#' y   <- dat$y
#' trt <- dat$trt
#'
#' ebal <- energy_balance(trt, x)
#'
#' print(ebal)
#'
#' # distribution of response:
#' quantile(y)
#'
#' # true trt effect:
#' dat$trt.eff
#'
#' # naive estimate of trt effect:
#' ipw_est(y, trt, rep(1, length(trt)))
#'
#' # estimated trt effect:
#' ipw_est(y, trt, ebal$weights)
#'
#' # estimated trt effect with true propensity:
#' wts_true <- 1 / (trt * dat$prob.trt + (1 - trt) * (1 - dat$prob.trt))
#' ipw_est(y, trt, wts_true)
#'
#' @export
energy_balance_2 <- function(trt,
                             x,
                             method      = c("ATE.3", "ATE", "ATT", "overlap"),
                             standardize = TRUE,
                             lambda      = 0,
                             verbose     = FALSE,
                             max.constr  = FALSE)
{

    type <- match.arg(method)

    if (standardize)
    {
        x.orig <- x
        x <- scale(x)
    }

    lambda <- lambda[1]
    stopifnot(lambda >= 0)

    #trt <- as.factor(as.vector(trt))

    if (!is.factor(trt))
    {
        trt <- as.factor(trt)
    }

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

    constr.sum  <- max.constr
    alpha       <- NULL
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
                trt_ind <- drop(1 * (trt == trt.levels[k]))
                QQ <- -trt_ind * t( trt_ind * t(QQ_all)) / n.vec[k] ^ 2

                aa <- 2 * as.vector(rowSums(trt_ind * QQ_all)) / (n.vec[k] * nn)
            } else
            {
                trt_ind <- drop(1 * (trt == trt.levels[k]))
                QQ <- QQ - trt_ind * t( trt_ind * t(QQ_all)) / n.vec[k] ^ 2

                aa <- aa + 2 * as.vector(rowSums(trt_ind * QQ_all)) / (n.vec[k] * nn)
            }
        }

        if (lambda > 0)
        {
            QQ <- QQ + lambda * diag(ncol(QQ)) / nn
        }


        rownames(QQ) <- paste(1:NROW(QQ))

        qf <- quadfun(Q = QQ, a = aa, id = rownames(QQ)) #quadratic obj.
        lb <- lbcon(0, id = rownames(QQ)) #lower bound

        AA           <- do.call(rbind, AA.list)
        rownames(AA) <- paste0("eq", 1:nrow(AA))
        sum.constr   <- n.vec

        AA_0 <- diag(nrow(QQ))

        Amat <- t(rbind(AA, AA_0))
        bvec <- c(n.vec, rep(0, nrow(QQ)))


        lc <- lincon(A   = AA,
                     dir = rep("==", length(sum.constr)),
                     val = sum.constr,
                     id  = rownames(QQ)) #linear constraint


        if (TRUE) # (lambda == 0)
        {
            if (constr.sum)
            {
                ub <- ubcon(10 * nn ^ (1/3), id = rownames(QQ))
                lcqp <- cop( f = qf, lb = lb, lc = lc, ub = ub)
            } else
            {
                lcqp <- cop( f = qf, lb = lb, lc = lc)
            }
        } else
        {
            if (constr.sum)
            {
                ub <- ubcon(10 * nn ^ (1/3), id = rownames(QQ))
                lcqp <- cop( f = qf, lb = lb, ub = ub)
            } else
            {
                lcqp <- cop( f = qf, lb = lb)
            }
        }

        # res <- solvecop(lcqp, solver = solver, quiet = !verbose)
        # print(isSymmetric(QQ))
        #
        # print(QQ[1:5,1:5])
        #
        #
        # QQold <- QQ
        #
        # QQ[upper.tri(QQ)] <- t(QQ)[upper.tri(QQ)]
        #
        # print(isSymmetric(QQ))
        #
        # print(QQ[1:5,1:5])
        #
        # QQ <- as.matrix(Matrix::forceSymmetric(QQ, uplo="L"))
        #
        # print(isSymmetric(QQ))
        #
        # #diag(QQ) <- 0.1
        #
        # Deig <- eigen(QQ, symmetric = TRUE)
        #
        # print(Deig$values)

        # res <- dykstra2(Dmat = QQ, dvec = -aa/2, Amat = Amat, bvec = bvec,
        #                 meq = 2, factorized = FALSE, eps = 1e-18)

        fn <- function(x)
        {
            val <- t(x) %*% QQ %*% x + t(aa) %*% x
            drop(val)
        }

        # res <- solnl(X = rep(1, nrow(QQ)),
        #              objfun = fn, Aeq = AA, Beq = n.vec,
        #              lb = rep(0, nrow(QQ)), ub = rep(nrow(QQ), nrow(QQ)))

        constrfun <- function(x)
        {
            drop(AA %*% x) - as.numeric(n.vec)
        }

        res <- slsqp(x0 = rep(1, nrow(QQ)), fn = fn,
                     lower = rep(0.0, nrow(QQ)),
                     heq = constrfun,
                     upper = rep(as.double(nrow(QQ)), nrow(QQ)))

        print(str(res))

        wts_quadopt <- drop(unname(res$par))
        print(wts_quadopt)

        print(sum(wts_quadopt[trt==1]))
        print(sum(wts_quadopt[trt!=1]))

        ### the overall objective function value
        final_value <- energy.dist.1bal.trt(wts_quadopt, x, trt, gamma = -1)

        ### the unweighted objective function value
        unweighted_value <- energy.dist.1bal.trt(rep(1, length(wts_quadopt)), x, trt, gamma = -1)

    } else if (type == "ATE.3")
    {
        ## three_way balances the treatment arm to the full population
        ##                    the control   arm to the full population
        ## and                the treatment arm to the control arm


        for (k in 1:K)
        {
            if (k == 1)
            {
                trt_ind <- drop(1 * (trt == trt.levels[k]))
                QQ <- -trt_ind * t( trt_ind * t(QQ_all)) / n.vec[k] ^ 2

                aa <- 2 * as.vector(rowSums(trt_ind * QQ_all)) / (n.vec[k] * nn)
            } else
            {
                trt_ind <- drop(1 * (trt == trt.levels[k]))
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

        if (lambda > 0)
        {
            QQ <- QQ + lambda * diag(ncol(QQ)) / nn
        }

        qf <- quadfun(Q = QQ, a = aa, id = rownames(QQ)) #quadratic obj.
        lb <- lbcon(0, id = rownames(QQ)) #lower bound


        lc <- lincon(A   = AA,
                     dir = rep("==", length(sum.constr)),
                     val = sum.constr,
                     id  = rownames(QQ)) #linear constraint


        if (TRUE) # (lambda == 0)
        {
            if (constr.sum)
            {
                ub <- ubcon(10 * nn ^ (1/3), id = rownames(QQ))
                lcqp <- cop( f = qf, lb = lb, lc = lc, ub = ub)
            } else
            {
                lcqp <- cop( f = qf, lb = lb, lc = lc)
            }
        } else
        {
            if (constr.sum)
            {
                ub <- ubcon(10 * nn ^ (1/3), id = rownames(QQ))
                lcqp <- cop( f = qf, lb = lb, ub = ub)
            } else
            {
                lcqp <- cop( f = qf, lb = lb)
            }
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

        if (lambda > 0)
        {
            QQ <- QQ + lambda * diag(ncol(QQ)) / nn
        }

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
        final_value <- weighted.energy.dist.trt2ctrl(unname(wts_quadopt), x, trt, gamma = -1)

        ### the unweighted objective function value
        unweighted_value <- weighted.energy.dist.trt2ctrl(rep(1, length(wts_quadopt)), x, trt, gamma = -1)
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

    ret
}



dykstra2 <- function (Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE,
          maxit = NULL, eps = NULL)
{
    Dmat <- as.matrix(Dmat)
    dvec <- as.numeric(dvec)
    pts <- length(dvec)
    if (nrow(Dmat) != pts | ncol(Dmat) != pts)
        stop("Inputs 'Dmat' and 'dvec' are incompatible.")
    Amat <- as.matrix(Amat)
    if (nrow(Amat) != pts)
        stop("Input 'Amat' must satisfy:  nrow(Amat) == length(dvec)")
    ncon <- ncol(Amat)
    if (missing(bvec)) {
        bvec <- rep(0, ncon)
    }
    else {
        bvec <- as.numeric(bvec)
        if (length(bvec) != ncon)
            stop("Input 'bvec' must satisfy:  length(bvec) == ncol(Amat)")
    }
    meq <- as.integer(meq[1])
    if (meq < 0 | meq > ncon)
        stop("Input 'meq' must be between 0 and length(bvec).")
    if (is.null(maxit)) {
        maxit <- 30 * pts
    }
    else {
        maxit <- as.integer(maxit[1])
        if (maxit < 1)
            stop("Input 'maxit' must be a positive integer.")
    }
    if (is.null(eps)) {
        eps <- .Machine$double.eps * pts
    }
    else {
        eps <- as.numeric(eps[1])
        if (eps < 0)
            stop("Input 'eps' must be a non-negative scalar.")
    }
    tol <- .Machine$double.eps * pts
    factorized <- as.logical(factorized[1])
    if (!factorized && !isSymmetric(Dmat)) {
        stop("Input 'Dmat' must be a symmetric matrix when 'factorized = FALSE'.")
    }
    uti <- upper.tri(Dmat)
    if (max(abs(Dmat[uti])) <= tol) {
        diag.flag <- TRUE
        if (factorized) {
            Rinv <- diag(Dmat)
        }
        else {
            Ddiag <- diag(Dmat)
            if (any(Ddiag < -tol))
                stop("Input 'Dmat' must be positive definite (or semidefinite).")
            nvals <- sum(Ddiag > (tol * max(Ddiag)))
            if (nvals < pts)
                Ddiag <- Ddiag + (tol * max(Ddiag) - min(Ddiag))
            Rinv <- 1/sqrt(Ddiag)
        }
    }
    else {
        diag.flag <- FALSE
        if (factorized) {
            Rinv <- Dmat
        }
        else {
            Deig <- eigen(Dmat, symmetric = TRUE)
            #if (any(Deig$values < -tol))
                #stop("Input 'Dmat' must be positive definite (or semidefinite).")
            nvals <- sum(Deig$values > (tol * Deig$values[1]))
            if (nvals < pts)
                Deig$values <- Deig$values + (tol * Deig$values[1] -
                                                  Deig$values[pts])
            Rinv <- matrix(0, pts, pts)
            for (i in 1:pts) Rinv[, i] <- Deig$vectors[, i]/sqrt(Deig$values[i])
        }
    }
    if (diag.flag) {
        gvec <- rep(0, pts)
        for (i in 1:pts) {
            gvec[i] <- Rinv[i] * dvec[i]
            Amat[i, ] <- Rinv[i] * Amat[i, ]
        }
    }
    else {
        gvec <- as.numeric(crossprod(Rinv, dvec))
        Amat <- crossprod(Rinv, Amat)
    }
    Acss <- colSums(Amat^2)
    if (diag.flag) {
        beta.unconstrained <- as.numeric(dvec/Rinv^2)
    }
    else {
        beta.unconstrained <- as.numeric(Rinv %*% crossprod(Rinv,
                                                            dvec))
    }
    maxb0 <- max(abs(beta.unconstrained))
    if (maxb0 > 1)
        eps <- eps * maxb0
    zeros <- rep(0, pts)
    beta.change <- matrix(0, pts, ncon)
    beta.solution <- beta.old <- gvec
    iter <- 0L
    ctol <- eps + 1
    while (ctol > eps && iter < maxit) {
        for (i in 1:ncon) {
            beta.work <- beta.solution - beta.change[, i]
            Ai <- sum(Amat[, i] * beta.work)
            passive <- ifelse(i <= meq, Ai == bvec[i], Ai >=
                                  bvec[i])
            if (passive) {
                beta.change[, i] <- zeros
                beta.solution <- beta.work
            }
            else {
                beta.change[, i] <- (bvec[i] - Ai) * Amat[, i]/Acss[i]
                beta.solution <- beta.work + beta.change[, i]
            }
        }
        ctol <- max(abs(beta.solution - beta.old))
        beta.old <- beta.solution
        iter <- iter + 1L
    }
    if (diag.flag) {
        beta.solution <- Rinv * beta.solution
    }
    else {
        beta.solution <- as.numeric(Rinv %*% beta.solution)
    }
    converged <- ifelse(ctol <= eps, TRUE, FALSE)
    if (factorized) {
        value <- NA
    }
    else {
        value.Q <- 0.5 * crossprod(beta.solution, Dmat %*% beta.solution)
        value.P <- sum(beta.solution * dvec)
        value <- as.numeric(value.Q - value.P)
    }
    results <- list(solution = beta.solution, value = value,
                    unconstrained = beta.unconstrained, iterations = iter,
                    converged = converged)
    class(results) <- "dykstra"
    return(results)
}
