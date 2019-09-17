

#' Simulation of confounded data
#'
#' @description Simulates confounded data
#'
#' @param n.obs number of observations
#' @param n.vars number of covariates
#' @param n.obs.test number of observations for test dataset
#' @param AR.cor correlation parameter for AR 1 correlation structure of covariates
#' @param y.model specification of model for response. Either \code{"A"}, \code{"B"}, or \code{"C"}
#' @param propensity.model specification of model for propensity score. Either \code{"I"}, \code{"II"}, or \code{"III"}
#'
#' @examples
#'
#' n <- 100
#' p <- 10
#'
#' set.seed(1)
#'
#' dat <- sim_confounded_data(n.obs = n, n.vars = p, AR.cor = 0.5,
#'                            propensity.model = "I", y.model = "I")
#'
#' x   <- dat$x
#' y   <- dat$y
#' trt <- dat$trt
#'
#' ncols <- 10
#' Cols <- rbPal(ncols)[as.numeric(cut(dat$prob.trt.test,breaks = ncols))]
#' plot(dat$x[,1:2], pch=19,cex=1,col=Cols,lwd=2,xlab="",ylab="",
#'      main = "propensity score vs X1 and X2")
#' legend("topleft",title="Decile",legend=levels(cut(dat$prob.trt.test,breaks = ncols)),col = rbPal(ncols),pch=19)
#'
#' plot(dat$x[,2:3], pch=19,cex=1,col=Cols,lwd=2,xlab="",ylab="",
#' main = "propensity score vs X1 and X2")
#'
#' plot(dat$x[,3:4], pch=19,cex=1,col=Cols,lwd=2,xlab="",ylab="",
#' main = "propensity score vs X1 and X2")
#' legend("bottomright",title="Decile",legend=levels(cut(dat$prob.trt.test,breaks = ncols)),col = rbPal(ncols),pch=19)
#'
#' plot(dat$x[,4:1], pch=19,cex=1,col=Cols,lwd=2,xlab="",ylab="",
#' main = "propensity score vs X1 and X2")
#' legend("bottomright",title="Decile",legend=levels(cut(dat$prob.trt.test,breaks = ncols)),col = rbPal(ncols),pch=19)
#'
#' @export
sim_confounded_data <- function(n.obs, n.vars = 10,
                                n.obs.test       = 10000,
                                AR.cor           = 0,
                                propensity.model = c("I", "II", "III", "IV", "V", "VI"),
                                y.model          = c("A", "B", "C", "D", "E", "F"),
                                transform.x      = FALSE)
{
    y.model          <- match.arg(y.model)
    propensity.model <- match.arg(propensity.model)

    if (n.vars < 5)
    {
        stop("'n.vars' must be 5 or greater")
    }

    stopifnot(abs(AR.cor) < 1)

    Sig      <- AR.cor ^ abs(outer(1:n.vars, 1:n.vars, FUN = "-"))
    x        <- mvrnorm(n = n.obs,      mu = rep(0, n.vars), Sigma = Sig)
    x.test   <- mvrnorm(n = n.obs.test, mu = rep(0, n.vars), Sigma = Sig)

    x.gen      <- x
    x.gen.test <- x.test

    if (transform.x)
    {
        x      <- transform.x.ks(x)
        x.test <- transform.x.ks(x.test)
    }

    y.model.func <- switch(y.model,
                           "A" = y.model.A,
                           "B" = y.model.B,
                           "C" = y.model.C,
                           "D" = y.model.D,
                           "E" = y.model.E,
                           "F" = y.model.F)

    prop.model.func <- switch(propensity.model,
                              "I"   = prop.model.I,
                              "II"  = prop.model.II,
                              "III" = prop.model.III,
                              "IV"  = prop.model.IV,
                              "V"   = prop.model.V,
                              "VI"  = prop.model.VI)

    prob_trt      <- prop.model.func(x.gen)
    prob_trt_test <- prop.model.func(x.gen.test)

    trt           <- rbinom(n.obs,      1, prob_trt)
    trt.test      <- rbinom(n.obs.test, 1, prob_trt_test)

    y_and_trt_effect_info      <- y.model.func(x.gen, trt)
    y_and_trt_effect_info.test <- y.model.func(x.gen.test, trt.test)

    y                          <- y_and_trt_effect_info$y
    trt_eff_sample             <- y_and_trt_effect_info$trt_eff_sample
    individual_trt_eff         <- y_and_trt_effect_info$individual_trt_eff

    y.test                     <- y_and_trt_effect_info.test$y
    trt_eff_sample.test        <- y_and_trt_effect_info.test$trt_eff_sample
    individual_trt_eff.test    <- y_and_trt_effect_info.test$individual_trt_eff

    if (!is.null(y_and_trt_effect_info$trt_eff))
    {
        trt_eff <- y_and_trt_effect_info$trt_eff
    } else
    {
        trt_eff <- mean(individual_trt_eff.test)
    }

    list(y = y, trt = trt, x = x, prob.trt = prob_trt,
         trt.eff = trt_eff,
         trt.eff.sample = trt_eff_sample, individual.trt.eff = individual_trt_eff,
         y.test = y.test, trt.test = trt.test, x.test = x.test, prob.trt.test = prob_trt_test,
         trt.eff.test = trt_eff_sample.test, individual.trt.eff.test = individual_trt_eff.test)
}


transform.x.ks <- function(x)
{
    X <- x

    X[,1] <- exp(x[,1] / 2)

    X[,2] <- x[,2] / (1 + exp(x[,1]) ) + 10

    X[,3] <- (x[,1] * x[,3] / 25 + 0.6 ) ^ 3

    X[,4] <- (20 + (x[,2] + x[,4]) ^ 2)

    X
}


prop.model.I <- function(x)
{
    ## orthogonalize first two variables
    x2 <- x; x2[,2] <- x[,2] - ((sum(x[,1] * x[,2])) / (sum(x[,1] ^ 2))) * x[,1]

    ## orthogonalize 3 and 4
    x2[,4] <- x[,4] - ((sum(x[,3] * x[,4])) / (sum(x[,3] ^ 2))) * x[,3]

    xbeta <- -0.5 + x2[,1] * x2[,2] - 0.5 * (x2[,1] - x2[,2]) ^ 2 + 1 * x2[,1] ^ 2 + 1 * x2[,1] ^ 2 -
        (x2[,3] * x2[,4] - 0.5 * (x2[,3] + x2[,4]) ^ 2 + 1 * x2[,3] ^ 2 + 1 * x2[,4] ^ 2)

    prob_trt <- 1 / (1 + exp(-xbeta))

    prob_trt
}


prop.model.II <- function(x)
{
    ## orthogonalize first two variables
    x2 <- x; x2[,2] <- x[,2] - ((sum(x[,1] * x[,2])) / (sum(x[,1] ^ 2))) * x[,1]

    ## orthogonalize 3 and 4
    x2[,4] <- x[,4] - ((sum(x[,3] * x[,4])) / (sum(x[,3] ^ 2))) * x[,3]


    xbeta <- -0.5 + x2[,5] * (x2[,1] * x2[,2] > 0 - 0.5 * (x2[,1] - x2[,2] > 0)) * (x2[,3] * x2[,4] > 0 - 0.5 * ((x2[,3] + x2[,4]) > 0) )

    prob_trt <- 1 / (1 + exp(-xbeta))

    prob_trt
}


prop.model.III <- function(x)
{

    x1x2f <- (x[,1] * x[,2]) ^ 2 > 0.5
    x1f   <- (x[,1])  ^ 2 > 0.1
    x2f   <- (x[,2])  ^ 2 > 0.1

    x3x4f <- (x[,3] * x[,4]) ^ 2 > 0.5
    x3f   <- (x[,3])  ^ 2 > 0.1
    x4f   <- (x[,4])  ^ 2 > 0.1

    xbeta <- -4 + 1 * (-2 * x1x2f + x1f + x2f) + 1 * (2 * x3x4f + x3f + x4f)

    prob_trt <- 1 / (1 + exp(-xbeta))

    prob_trt
}

prop.model.IV <- function(x)
{

    xbeta <- numeric(NROW(x))

    nv <- 4
    ct <- 0
    for (j in 1:(nv - 1))
    {
        for (k in (j:nv))
        {
            ct <- ct + 1

            if (ct %% 3 == 0)
            {
                sgn <- 1
            } else
            {
                sgn <- -1
            }

            xbeta <- xbeta + sgn * (x[,j] * x[,k])
        }
    }

    xbeta <- 10 * xbeta / var(xbeta)

    prob_trt <- 1 / (1 + exp(-xbeta))

    prob_trt
}



prop.model.V <- function(x)
{
    ## cross-in-tray function
    ## https://www.sfu.ca/~ssurjano/drop.html

    x1 <- x[,1]
    x2 <- x[,2]
    x3 <- x[,3]
    x4 <- x[,4]

    x12_diff <- x1 - x2
    x23_diff <- x2 - x3
    x34_diff <- x3 - x4


    xbeta <- x12_diff ^ 2 - x23_diff ^ 2 + x34_diff * x1 * x2

    prob_trt <- 1 / (1 + exp(-xbeta))

    prob_trt
}





prop.model.VI <- function(x)
{

    xbeta <- 0.5 * (x[,1] - x[,2] + x[,3] - x[,4] + x[,5])

    xbeta <- -x[,1] + 0.5 * x[,2] - 0.25 * x[,3] - 0.1 * x[,4]

    prob_trt <- 1 / (1 + exp(-xbeta))

    prob_trt
}

y.model.A <- function(x, trt, sd.error = 1)
{
    trt_eff <- NULL

    E_y_given_x_trt <- function(xx, trtt)
    {

        #E_y <- 210 + 10 * trtt + 2 * xx[,1] ^ 2 + 2 * xx[,2]^2 - 1 / (1 + (xx[,1] * xx[,2]) ^ 2)  - 5 * xx[,3] * xx[,4] - 5 * xx[,1] * xx[,2]

        ## friedman function
        xx <- apply(xx, 2, function(xc) (xc - min(xc)) / (max(xc) - min(xc)) )
        E_y <- (trtt + 1) * (10 * sin(pi * xx[,1] * xx[,2]) + 20 * (xx[,3] - 0.5) ^ 2 + 10 * xx[,4] + 5 * xx[,5])
        E_y
    }

    n <- NROW(x)

    y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)

    individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)

    trt_eff_sample <- mean(individual_trt_eff)

    list(y = y,
         trt_eff_sample     = trt_eff_sample,
         individual_trt_eff = individual_trt_eff,
         trt_eff = trt_eff)
}

y.model.B <- function(x, trt, sd.error = 1)
{
    trt_eff <- NULL

    E_y_given_x_trt <- function(xx, trtt)
    {
        E_y <- 210 + (1 * trtt + 0.5) * (27.4 * xx[,1] + 13.7 * xx[,2] + 13.7 * xx[,3] + 13.7 * xx[,4])
        E_y
    }

    n <- NROW(x)

    y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)

    individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)

    trt_eff_sample <- mean(individual_trt_eff)

    list(y = y,
         trt_eff_sample     = trt_eff_sample,
         individual_trt_eff = individual_trt_eff,
         trt_eff = trt_eff)
}

y.model.C <- function(x, trt, sd.error = 1)
{
    trt_eff <- NULL

    E_y_given_x_trt <- function(xx, trtt)
    {
        E_y <- 210 + (2 * trtt + 1) * (5 * xx[,1] * xx[,2] + 5 * xx[,3] * xx[,4] - 5 * (xx[,1] - xx[,2]) ^ 2 + 5 * (xx[,3] - xx[,4]) ^ 2  )
        E_y
    }

    n <- NROW(x)

    y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)

    individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)

    trt_eff_sample <- mean(individual_trt_eff)

    list(y = y,
         trt_eff_sample     = trt_eff_sample,
         individual_trt_eff = individual_trt_eff,
         trt_eff = trt_eff)
}

y.model.D <- function(x, trt, sd.error = 1)
{
    trt_eff <- NULL

    E_y_given_x_trt <- function(xx, trtt)
    {
        E_y <- 210 + (27.4 * xx[,1] + 13.7 * xx[,2] + 13.7 * xx[,3] + 13.7 * xx[,4])
        E_y
    }

    n <- NROW(x)

    y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)

    individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)

    trt_eff <- 0

    trt_eff_sample <- mean(individual_trt_eff)

    list(y = y,
         trt_eff_sample     = trt_eff_sample,
         individual_trt_eff = individual_trt_eff,
         trt_eff = trt_eff)
}


y.model.E <- function(x, trt, sd.error = 1)
{
    trt_eff <- NULL

    nv <- 4

    E_y_given_x_trt <- function(xx, trtt)
    {

        xbeta <- numeric(NROW(x))

        ct <- 0
        for (j in 1:(nv - 1))
        {
            for (k in (j:nv))
            {
                ct <- ct + 1

                if (ct %% 2 == 0)
                {
                    sgn <- 1
                } else
                {
                    sgn <- -1
                }

                xbeta <- xbeta + sgn * (xx[,j] * xx[,k] - 0.25 * (xx[,j] + xx[,k]))
            }
        }

        E_y <- 10 * (trtt + 1) * xbeta / var(xbeta)

        E_y
    }

    n <- NROW(x)

    y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)

    individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)

    trt_eff_sample <- mean(individual_trt_eff)

    list(y = y,
         trt_eff_sample     = trt_eff_sample,
         individual_trt_eff = individual_trt_eff,
         trt_eff = trt_eff)
}


y.model.F <- function(x, trt, sd.error = 1)
{
    trt_eff <- NULL

    E_y_given_x_trt <- function(xx, trtt)
    {

        x1 <- xx[,1]
        x2 <- xx[,2]
        x3 <- xx[,3]
        x4 <- xx[,4]

        x12_diff <- x1 - x2
        x23_diff <- x2 - x3
        x34_diff <- x3 - x4
        x41_diff <- x4 - x1


        E_y <- 2 * (trtt + 0.5) * exp(-x12_diff) * sin(x3 ^ 2) * cos(x4 ^ 2) +
            2 * (trtt + 0.5) * exp(-x23_diff) * sin(x1 ^ 2) * cos(x4 ^ 2) +
            exp(-x34_diff) * (x3) * x4 + exp(-x41_diff) * (x2) * x3

        E_y
    }

    n <- NROW(x)

    y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)

    individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)


    trt_eff_sample <- mean(individual_trt_eff)

    list(y = y,
         trt_eff_sample     = trt_eff_sample,
         individual_trt_eff = individual_trt_eff,
         trt_eff = trt_eff)
}



y.model.G <- function(x, trt, sd.error = 1)
{
    trt_eff <- NULL

    E_y_given_x_trt <- function(xx, trtt)
    {
        #E_y <- 210 + (0.5 * trtt + 0.5) * (2 * xx[,1] ^ 2 + 2 * xx[,2]^2 - 1 / (1 + (xx[,1] * xx[,2]) ^ 2)  - 1 * (xx[,3] * xx[,4])^3 - 1 * (xx[,1] * xx[,2]) ^ 3)
        E_y <- 210 + (0.5 * trtt + 0.5) * (2 * xx[,1] ^ 2 + 2 * xx[,2]^2 - 1 * (xx[,3] * xx[,4])^3 - 1 * (xx[,1] * xx[,2]) ^ 3)
        E_y
    }

    n <- NROW(x)

    y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)

    individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)


    trt_eff_sample <- mean(individual_trt_eff)

    list(y = y,
         trt_eff_sample     = trt_eff_sample,
         individual_trt_eff = individual_trt_eff,
         trt_eff = trt_eff)
}




