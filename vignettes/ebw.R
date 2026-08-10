## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4)
library(ebw)
set.seed(20260809)

## ----binary-fit---------------------------------------------------------------
dat <- sim_confounded_data(n.obs = 300, n.vars = 6, AR.cor = 0.3)
fit_b <- energy_balance(x = dat$x, treatment = dat$trt, estimand = "ATE")
fit_b

## ----binary-weights-----------------------------------------------------------
w <- weights(fit_b)
str(w)

## ----binary-predict-----------------------------------------------------------
new <- sim_confounded_data(n.obs = 50, n.vars = 6, AR.cor = 0.3)
w_new <- predict(fit_b, newdata = new$x)
summary(w_new)

## ----binary-balance-----------------------------------------------------------
balance(fit_b)
balance(fit_b, newdata = list(x = new$x, trt = new$trt))

## ----binary-plot, fig.alt = "Love plot of covariate balance"------------------
plot(fit_b, type = "balance")
plot(fit_b, type = "weights")

## ----cobalt, eval = requireNamespace("cobalt", quietly = TRUE)----------------
cobalt::bal.tab(dat$x, treat = dat$trt, weights = weights(fit_b))

## ----cont-fit-----------------------------------------------------------------
cdat <- sim_continuous_data(n.obs = 300, n.vars = 4, alpha.norm2 = 0.5)
fit_c <- independence_weights(A = cdat$a, X = cdat$x)
fit_c

## ----cont-predict-------------------------------------------------------------
cnew <- sim_continuous_data(n.obs = 50, n.vars = 4, alpha.norm2 = 0.5)
w_cnew <- predict(fit_c, newdata = list(x = cnew$x, a = cnew$a))
summary(w_cnew)

## ----cont-balance-------------------------------------------------------------
balance(fit_c)

## ----cont-constrained---------------------------------------------------------
fit_pm <- independence_weights(A = cdat$a, X = cdat$x, preserve_means = TRUE)
w_pm <- weights(fit_pm)

# the weighted covariate means now match the unweighted means
round(max(abs(colSums(w_pm * cdat$x) / sum(w_pm) - colMeans(cdat$x))), 12)

# and prediction remains available out of sample
head(predict(fit_pm, newdata = list(x = cnew$x, a = cnew$a)))

