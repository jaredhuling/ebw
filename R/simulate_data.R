

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
                                propensity.model = c("I", "II", "III", "IV", "V", "VI", "VII"),
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

    # if (propensity.model == "VI")
    # {
    #     sd_vec  <- runif(n.obs, min = 0.25, 3)
    #     sd_vec2 <- runif(n.obs, min = 0.25, 3)
    #     cov_vec <- pmin(sd_vec, sd_vec2)
    #
    #     param_list <- list(sd_vec = sd_vec, sd_vec2 = sd_vec, cov_vec = cov_vec)
    #
    #
    #     sd_vec  <- runif(n.obs.test, min = 0.25, 3)
    #     sd_vec2 <- runif(n.obs.test, min = 0.25, 3)
    #     cov_vec <- pmin(sd_vec, sd_vec2)
    #
    #     param_test_list <- list(sd_vec = sd_vec, sd_vec2 = sd_vec, cov_vec = cov_vec)
    # } else
    # {
    #     param_list <- param_test_list <- NULL
    # }

    param_list <- param_test_list <- NULL

    y.model.func <- switch(y.model,
                           "A" = y.model.A,
                           "B" = y.model.B,
                           "C" = y.model.C,
                           "D" = y.model.D,
                           "E" = y.model.E,
                           "F" = y.model.F)

    prop.model.func <- switch(propensity.model,
                              "I"   = ebw:::prop.model.I,
                              "II"  = ebw:::prop.model.II,
                              "III" = ebw:::prop.model.III,
                              "IV"  = ebw:::prop.model.IV,
                              "V"   = ebw:::prop.model.V,
                              "VI"  = ebw:::prop.model.VI,
                              "VII" = ebw:::prop.model.VII)

    prob_trt      <- prop.model.func(x.gen, param_list)
    prob_trt_test <- prop.model.func(x.gen.test, param_test_list)

    trt           <- rbinom(n.obs,      1, prob_trt)
    trt.test      <- rbinom(n.obs.test, 1, prob_trt_test)

    y_and_trt_effect_info      <- y.model.func(x.gen, trt)
    y_and_trt_effect_info.test <- y.model.func(x.gen.test, trt.test)

    if (transform.x)
    {
        x      <- transform.x.ks(x)
        x.test <- transform.x.ks(x.test)
    }

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
         trt.eff.sample = trt_eff_sample,
         individual.trt.eff = individual_trt_eff,
         y.test = y.test, trt.test = trt.test,
         x.test = x.test, prob.trt.test = prob_trt_test,
         trt.eff.test = trt_eff_sample.test,
         individual.trt.eff.test = individual_trt_eff.test)
}


transform.x.ks <- function(x)
{
    X <- x

    X[,1] <- exp(x[,1] / 2)

    X[,2] <- x[,2] / (1 + exp(x[,1]) ) + 10

    X[,3] <- (x[,1] * x[,3] / 25 + 0.6 ) ^ 3

    X[,4] <- (20 + (x[,2] + x[,4]) ^ 2)

    if (ncol(X) >= 8)
    {
        X[,5] <- exp(-x[,5] / 2)

        X[,6] <- x[,6] / (1 + exp(x[,5]) ) + 10

        X[,7] <- (x[,1] * x[,7] / 25 + 0.6 ) ^ 3

        X[,8] <- (5 + (x[,6] + x[,8]) ^ 2)
    }

    X
}











prop.model.I <- function(x, param_list = NULL)
{

    x1 <- x[,1]
    x2 <- x[,2]
    x3 <- x[,3]
    x4 <- x[,4]

    x12_diff <- (x1 - x2)
    x23_diff <- (x2 - x3)
    x34_diff <- (x3 - x4)


    xbeta <- 2 * x1 * x2 * (abs(x1) > 1) * (abs(x2) > 1) +
        2 * x2 * x3 * (abs(x2) < 1) * (abs(x3) < 1) +
        2 * x3 * x4 * (abs(x3) > 1) * (abs(x4) > 1) +
        2 * x1 * x4 * (abs(x1) < 1) * (abs(x4) < 1) +
        1 * (abs(x1) > 0.5) * (abs(x2) > 0.5) * (abs(x3) > 0.5) * (abs(x4) > 0.5) +
        1 * (abs(x1) < 0.25) * (abs(x2) > 0.25) * (abs(x3) < 0.25) * (abs(x4) > 0.25)

    prob_trt <- 1 / (1 + exp(-xbeta))

    prob_trt
}


prop.model.II <- function(x, param_list = NULL)
{

    x1 <- x[,1]
    x2 <- x[,2]
    x3 <- x[,3]
    x4 <- x[,4]

    x12_diff <- x1 - x2
    x23_diff <- x2 - x3
    x34_diff <- x3 - x4


    xbeta <- -2 + log(abs(x12_diff)) - log(abs(x23_diff)) + abs(x34_diff * x1 * x2) ^ (1/2)

    prob_trt <- 1 / (1 + exp(-xbeta))

    prob_trt
}



prop.model.III <- function(x, param_list = NULL)
{

    xbeta <- -x[,1] + 0.5 * x[,2] - 0.25 * x[,3] - 0.1 * x[,4] -
              x[,5] + 0.5 * x[,6] - 0.25 * x[,7] - 0.1 * x[,8]

    prob_trt <- 1 / (1 + exp(-xbeta))

    prob_trt
}

prop.model.IV <- function(x, param_list = NULL)
{

    xbeta <- numeric(NROW(x))

    nv <- 4
    ct <- 0
    for (j in 1:(nv - 1))
    {
        for (k in (j:nv))
        {
            ct <- ct + 1

            # if (ct %% 3 == 0)
            # {
            #     sgn <- 1
            # } else
            # {
            #     sgn <- -1
            # }

            sgn <- (-1) ^ (2*k-j)

            xbeta <- xbeta + sgn * (x[,j] * x[,k])
        }
    }

    xbeta <- 5 * xbeta / var(xbeta)

    prob_trt <- 1 / (1 + exp(-xbeta))

    prob_trt
}


prop.model.V <- function(x, param_list = NULL)
{

    xbeta <- -2 + 2*(x[,1] * x[,2] + 0.5 * (x[,1] - x[,2]) ^ 2 -
        (x[,3] * x[,4] + 0.5 * (x[,3] + x[,4]) ^ 2))

    pi <- 1 / (1 + exp(-xbeta))

    pi
}




prop.model.VII <- function(x, param_list = NULL)
{
    x1 <- x[,1]
    x2 <- x[,2]
    x3 <- x[,3]
    x4 <- x[,4]
    x5 <- x[,5]
    x6 <- x[,6]
    x7 <- x[,7]
    x8 <- x[,8]

    x12_diff <- x1 - 2 * x2
    x23_diff <- x2 - 2 * x3
    x34_diff <- x3 - 2 * x4
    x45_diff <- x4 - 2 * x5

    xbeta <- abs(x12_diff) * abs(x23_diff) - abs(x34_diff) * abs(x45_diff) + x6 - 0.5 * x7 - 0.25 * x8

    prob_trt <- exp(xbeta) / (1 + exp(xbeta))

    prob_trt
}


prop.model.VI <- function(x, param_list = NULL)
{
    x1 <- x[,1]
    x2 <- x[,2]
    x3 <- x[,3]
    x4 <- x[,4]
    x5 <- x[,5]
    x6 <- x[,6]
    x7 <- x[,7]
    x8 <- x[,8]

    x12_diff <- x1 - 2 * x2
    x23_diff <- x2 - 2 * x3
    x34_diff <- x3 - 2 * x4
    x45_diff <- x4 - 2 * x5

    xbeta <- (x12_diff) * (x23_diff) - (x34_diff) * (x45_diff) + x6 - 0.5 * x7 - 0.25 * x8

    prob_trt <- exp(xbeta) / (1 + exp(xbeta))

    prob_trt
}


#### OUTCOME MODELS




y.model.A <- function(x, trt, sd.error = 1)
{
    trt_eff <- NULL

    E_y_given_x_trt <- function(xx, trtt)
    {
        E_y <- 210 + (27.4 * abs(xx[,1]) + 13.7 * abs(xx[,2]) + 13.7 * abs(xx[,3]) + 13.7 * abs(xx[,4]))
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



y.model.B <- function(x, trt, sd.error = 1)
{
    trt_eff <- NULL

    E_y_given_x_trt <- function(xx, trtt)
    {
        E_y <- xx[,1] * (xx[,2] ^ 3) * (xx[,3] ^ 2) * xx[,4] + xx[,4] * abs(xx[,1]) ^ (1/2)
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

        x1 <- xx[,1]
        x2 <- xx[,2]
        x3 <- xx[,3]
        x4 <- xx[,4]
        x5 <- xx[,5]
        x6 <- xx[,6]
        x7 <- xx[,7]
        x8 <- xx[,8]
        x9 <- xx[,9]
        x10 <- xx[,10]

        x12_diff <- x1 - 2 * x2
        x23_diff <- x2 - 2 * x3
        x34_diff <- x3 - 2 * x4
        x45_diff <- x4 - 2 * x5


        E_y <- c(2 - 2 * x1*(x1>0) * trtt) * (x12_diff) +
            c(2 - 2 * x2*(x2>0) * trtt) * (x23_diff) +
            c(2 - 2 * x3*(x3>0) * trtt) * (x34_diff) +
            c(2 - 2 * x4*(x4>0) * trtt) * (x45_diff)

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

        x1 <- xx[,1]
        x2 <- xx[,2]
        x3 <- xx[,3]
        x4 <- xx[,4]
        x5 <- xx[,5]
        x6 <- xx[,6]
        x7 <- xx[,7]
        x8 <- xx[,8]
        x9 <- xx[,9]
        x10 <- xx[,10]

        beta <- c(0.8, 0.25, 0.6, -0.4, -0.8, -0.5, 0.7)

        f_A <- drop(xx[,1:7] %*% beta)

        f_G <- f_A + beta[2] * x2^2 + beta[4] * x4^2 + beta[7] * x7^2 + 0.5 * beta[1]*x1*x3 + 0.7*beta[2]*x2*x4 +
            0.5*beta[3]*x3*x5 + 0.7*beta[4]*x4*x6 + 0.5*beta[5]*x5*x7 + 0.5*beta[1]*x1*x6+0.7*beta[2]*x2*x3+0.5*beta[3]*x3*x4 +
            0.5*beta[4]*x4*x5 + 0.5*beta[5]*x5*x6

        E_y <- 5 * f_G

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



y.model.E <- function(x, trt, sd.error = 1)
{
    trt_eff <- NULL

    E_y_given_x_trt <- function(xx, trtt)
    {
        E_y <- 210 + (1.5 * trtt - 0.5) * (27.4 * xx[,1] + 13.7 * xx[,2] + 13.7 * xx[,3] + 13.7 * xx[,4])
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


# prop.model.I <- function(x, param_list = NULL)
# {
#     ## orthogonalize first two variables
#     x2 <- x; x2[,2] <- x[,2] - ((sum(x[,1] * x[,2])) / (sum(x[,1] ^ 2))) * x[,1]
#
#     ## orthogonalize 3 and 4
#     x2[,4] <- x[,4] - ((sum(x[,3] * x[,4])) / (sum(x[,3] ^ 2))) * x[,3]
#
#     xbeta <- -0.5 + x2[,1] * x2[,2] - 0.5 * (x2[,1] - x2[,2]) ^ 2 + 1 * x2[,1] ^ 2 + 1 * x2[,1] ^ 2 -
#         (x2[,3] * x2[,4] - 0.5 * (x2[,3] + x2[,4]) ^ 2 + 1 * x2[,3] ^ 2 + 1 * x2[,4] ^ 2)
#
#     prob_trt <- 1 / (1 + exp(-xbeta))
#
#     prob_trt
# }
#
#
# prop.model.II <- function(x, param_list = NULL)
# {
#     ## orthogonalize first two variables
#     x2 <- x; x2[,2] <- x[,2] - ((sum(x[,1] * x[,2])) / (sum(x[,1] ^ 2))) * x[,1]
#
#     ## orthogonalize 3 and 4
#     x2[,4] <- x[,4] - ((sum(x[,3] * x[,4])) / (sum(x[,3] ^ 2))) * x[,3]
#
#
#     xbeta <- -0.5 + x2[,5] * ( (x2[,1] * x2[,2] > 0) - 0.5 * (x2[,1] - x2[,2] > 0)) * ( (x2[,3] * x2[,4] > 0) - 0.5 * ((x2[,3] + x2[,4]) > 0) )
#
#     prob_trt <- 1 / (1 + exp(-xbeta))
#
#     prob_trt
# }
#
#
# prop.model.III <- function(x, param_list = NULL)
# {
#
#     x1x2f <- (x[,1] * x[,2]) ^ 2 > 0.5
#     x1f   <- (x[,1])  ^ 2 > 0.1
#     x2f   <- (x[,2])  ^ 2 > 0.1
#
#     x3x4f <- (x[,3] * x[,4]) ^ 2 > 0.5
#     x3f   <- (x[,3])  ^ 2 > 0.1
#     x4f   <- (x[,4])  ^ 2 > 0.1
#
#     xbeta <- -4 + 1 * (-2 * x1x2f + x1f + x2f) + 1 * (2 * x3x4f + x3f + x4f)
#
#     prob_trt <- 1 / (1 + exp(-xbeta))
#
#     prob_trt
# }
#
# prop.model.IV <- function(x, param_list = NULL)
# {
#
#     xbeta <- numeric(NROW(x))
#
#     nv <- 4
#     ct <- 0
#     for (j in 1:(nv - 1))
#     {
#         for (k in (j:nv))
#         {
#             ct <- ct + 1
#
#             if (ct %% 3 == 0)
#             {
#                 sgn <- 1
#             } else
#             {
#                 sgn <- -1
#             }
#
#             xbeta <- xbeta + sgn * (x[,j] * x[,k])
#         }
#     }
#
#     xbeta <- 10 * xbeta / var(xbeta)
#
#     prob_trt <- 1 / (1 + exp(-xbeta))
#
#     prob_trt
# }
#
#
#
# prop.model.V <- function(x, param_list = NULL)
# {
#
#     x1 <- x[,1]
#     x2 <- x[,2]
#     x3 <- x[,3]
#     x4 <- x[,4]
#
#     x12_diff <- x1 - x2
#     x23_diff <- x2 - x3
#     x34_diff <- x3 - x4
#
#
#     xbeta <- x12_diff ^ 2 - x23_diff ^ 2 + x34_diff * x1 * x2
#
#     prob_trt <- 1 / (1 + exp(-xbeta))
#
#     prob_trt
# }
#
#
# prop.model.VI <- function(x, param_list = NULL)
# {
#     x1 <- x[,1]
#     x2 <- x[,2]
#     x3 <- x[,3]
#     x4 <- x[,4]
#     x5 <- x[,5]
#     x6 <- x[,6]
#     x7 <- x[,7]
#     x8 <- x[,8]
#
#     x12_diff <- x1 - 2 * x2
#     x23_diff <- x2 - 2 * x3
#     x34_diff <- x3 - 2 * x4
#     x45_diff <- x4 - 2 * x5
#
#     xbeta <- 0 +
#         1 * ( (x1 > 2) + (x2 > 2) - (x3 > 2) - (x4 > 2) + (x5 > 1) - (x6 > 1)) - 0.5 * x5 ^ 2 + 0.5 * x6 ^ 2 +
#         0.5 * x6 - 0.25 * x7 - 0.1 * x8
#
#     prob_trt <- exp(xbeta) / (1 + exp(xbeta))
#
#     prob_trt
# }
#
#
# prop.model.VII <- function(x, param_list = NULL)
# {
#
#     xbeta <- -x[,1] + 0.5 * x[,2] - 0.25 * x[,3] - 0.1 * x[,4]
#
#     prob_trt <- 1 / (1 + exp(-xbeta))
#
#     prob_trt
# }
#
#
#
# y.model.A <- function(x, trt, sd.error = 1)
# {
#     trt_eff <- NULL
#
#     E_y_given_x_trt <- function(xx, trtt)
#     {
#         E_y <- xx[,1] * (xx[,2] ^ 3) * (xx[,3] ^ 2) * xx[,4] + xx[,1] * abs(xx[,1]) ^ (1/2)
#         E_y
#     }
#
#     n <- NROW(x)
#
#     y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)
#
#     individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)
#
#     trt_eff_sample <- mean(individual_trt_eff)
#
#     list(y = y,
#          trt_eff_sample     = trt_eff_sample,
#          individual_trt_eff = individual_trt_eff,
#          trt_eff = trt_eff)
# }
#
# y.model.B <- function(x, trt, sd.error = 1)
# {
#     trt_eff <- NULL
#
#     E_y_given_x_trt <- function(xx, trtt)
#     {
#         E_y <- 210 + (1 * trtt + 0.5) * (27.4 * xx[,1] + 13.7 * xx[,2] + 13.7 * xx[,3] + 13.7 * xx[,4])
#         E_y
#     }
#
#     n <- NROW(x)
#
#     y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)
#
#     individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)
#
#     trt_eff_sample <- mean(individual_trt_eff)
#
#     list(y = y,
#          trt_eff_sample     = trt_eff_sample,
#          individual_trt_eff = individual_trt_eff,
#          trt_eff = trt_eff)
# }
#
# y.model.C <- function(x, trt, sd.error = 1)
# {
#     trt_eff <- NULL
#
#     E_y_given_x_trt <- function(xx, trtt)
#     {
#         E_y <- 210 + (2 * trtt + 1) * (5 * xx[,1] * xx[,2] + 5 * xx[,3] * xx[,4] - 5 * (xx[,1] - xx[,2]) ^ 2 + 5 * (xx[,3] - xx[,4]) ^ 2  )
#         E_y
#     }
#
#     n <- NROW(x)
#
#     y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)
#
#     individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)
#
#     trt_eff_sample <- mean(individual_trt_eff)
#
#     list(y = y,
#          trt_eff_sample     = trt_eff_sample,
#          individual_trt_eff = individual_trt_eff,
#          trt_eff = trt_eff)
# }
#
# y.model.D <- function(x, trt, sd.error = 1)
# {
#     trt_eff <- NULL
#
#     E_y_given_x_trt <- function(xx, trtt)
#     {
#         E_y <- 210 + (27.4 * xx[,1] + 13.7 * xx[,2] + 13.7 * xx[,3] + 13.7 * xx[,4])
#         E_y
#     }
#
#     n <- NROW(x)
#
#     y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)
#
#     individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)
#
#     trt_eff <- 0
#
#     trt_eff_sample <- mean(individual_trt_eff)
#
#     list(y = y,
#          trt_eff_sample     = trt_eff_sample,
#          individual_trt_eff = individual_trt_eff,
#          trt_eff = trt_eff)
# }
#
#
# y.model.E <- function(x, trt, sd.error = 1)
# {
#     trt_eff <- NULL
#
#     nv <- 4
#
#     E_y_given_x_trt <- function(xx, trtt)
#     {
#
#         xbeta <- numeric(NROW(x))
#
#         ct <- 0
#         for (j in 1:(nv - 1))
#         {
#             for (k in (j:nv))
#             {
#                 ct <- ct + 1
#
#                 if (ct %% 2 == 0)
#                 {
#                     sgn <- 1
#                 } else
#                 {
#                     sgn <- -1
#                 }
#
#                 xbeta <- xbeta + sgn * (xx[,j] * xx[,k] - 0.25 * (xx[,j] + xx[,k]))
#             }
#         }
#
#         E_y <- 10 * (trtt + 1) * xbeta / var(xbeta)
#
#         E_y
#     }
#
#     n <- NROW(x)
#
#     y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)
#
#     individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)
#
#     trt_eff_sample <- mean(individual_trt_eff)
#
#     list(y = y,
#          trt_eff_sample     = trt_eff_sample,
#          individual_trt_eff = individual_trt_eff,
#          trt_eff = trt_eff)
# }
#
#
# y.model.F <- function(x, trt, sd.error = 1)
# {
#     trt_eff <- NULL
#
#     E_y_given_x_trt <- function(xx, trtt)
#     {
#
#         x1 <- xx[,1]
#         x2 <- xx[,2]
#         x3 <- xx[,3]
#         x4 <- xx[,4]
#         x5 <- xx[,5]
#         x6 <- xx[,6]
#         x7 <- xx[,7]
#         x8 <- xx[,8]
#         x9 <- xx[,9]
#         x10 <- xx[,10]
#
#         x12_diff <- x1 - 2 * x2
#         x23_diff <- x2 - 2 * x3
#         x34_diff <- x3 - 2 * x4
#         x45_diff <- x4 - 2 * x5
#
#
#         E_y <- 210 + 2.5 * trtt * (x12_diff) / (0.2 + (x3) ^ 2 ) +
#             2.5 * (1 - trtt) * (x23_diff) / (0.2 + (x4) ^ 2 ) +
#             10 * (x34_diff) / (0.2 + (x5) ^ 2 ) +
#             10 * (x45_diff) / (0.2 + (x6) ^ 2 ) +
#             10 * (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
#
#         E_y
#     }
#
#     n <- NROW(x)
#
#     y <- E_y_given_x_trt(x, trt) + rnorm(n, sd = sd.error)
#
#     individual_trt_eff <- E_y_given_x_trt(x, 1) - E_y_given_x_trt(x, 0)
#
#
#     trt_eff_sample <- mean(individual_trt_eff)
#
#     list(y = y,
#          trt_eff_sample     = trt_eff_sample,
#          individual_trt_eff = individual_trt_eff,
#          trt_eff = trt_eff)
# }
#
