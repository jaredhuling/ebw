# Grouped energy balancing: standard / improved binary and multi-category,
# out-of-sample prediction, backward compatibility, and the effect of the
# between-group (improved) term.

# Shared data.
make_data <- function() {
  set.seed(20260810L)
  n <- 300; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  tb <- rbinom(n, 1, plogis(X[, 1] - X[, 2]))            # binary 0/1
  lin <- X %*% c(1, -1, 0.5, 0, 0)
  pr <- exp(cbind(0, lin, -lin)); pr <- pr / rowSums(pr)
  tm <- apply(pr, 1, function(pp) sample(1:3, 1, prob = pp))  # 3-level
  list(X = X, tb = tb, tm = tm, n = n, p = p)
}

# Self-consistency helper: predict at the training points (each carrying its own
# treatment) must reproduce the fitted QP weights.
self_consistency <- function(fit) {
  pred <- predict(fit, newdata = fit$x, treatment = fit$treatment)
  max(abs(pred - fit$weights))
}

group_sums <- function(fit) {
  levs <- fit$treatment_levels
  vapply(levs, function(lv) sum(fit$weights[fit$treatment == lv]), numeric(1))
}

test_that("self-consistency holds for all four grouped variants", {
  d <- make_data()

  f_std <- energy_balance(x = d$X, treatment = d$tb)                    # binary beta=0
  f_imp <- energy_balance(x = d$X, treatment = d$tb, improved = TRUE)   # binary beta=1
  f_m0  <- energy_balance(x = d$X, treatment = d$tm)                    # K=3 beta=0
  f_m1  <- energy_balance(x = d$X, treatment = d$tm, improved = TRUE)   # K=3 beta=1

  expect_lt(self_consistency(f_std), 1e-4)
  expect_lt(self_consistency(f_imp), 1e-4)
  expect_lt(self_consistency(f_m0),  1e-4)
  expect_lt(self_consistency(f_m1),  1e-4)

  # each group's weights sum to one
  for (fit in list(f_std, f_imp, f_m0, f_m1))
    expect_true(all(abs(group_sums(fit) - 1) < 1e-6))

  expect_equal(f_std$treatment_type, "binary")
  expect_equal(f_imp$treatment_type, "binary")
  expect_equal(f_m0$treatment_type, "multinomial")
  expect_equal(f_m0$n_groups, 3L)
})

test_that("the improved term reduces the between-group energy distance", {
  d <- make_data()
  f_std <- energy_balance(x = d$X, treatment = d$tb)
  f_imp <- energy_balance(x = d$X, treatment = d$tb, improved = TRUE)

  btw_standard <- balance(f_std)$between_group$adjusted
  btw_improved <- balance(f_imp)$between_group$adjusted
  expect_length(btw_standard, 1L)
  expect_lt(btw_improved, btw_standard)
})

test_that("standard binary is unchanged from the pre-change engine", {
  set.seed(20260810L)
  n <- 250; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  tb <- rbinom(n, 1, plogis(X[, 1] - X[, 2]))

  set.seed(99L)
  fit <- energy_balance(x = X, treatment = tb)     # default: binary, beta = 0

  # replicate the internal binary path exactly (same scaled X, same RNG stream)
  sc <- .scale_covariates(X, "std")
  set.seed(99L)
  eng <- .fit_binary(sc$X, tb, estimand = "ATE", level = NULL, eps = 1e-2,
                     kernel = energy_kernel_X, max_iter = 20000L, tol = 1e-11)

  expect_lt(max(abs(fit$weights - eng$weights)), 1e-8)
})

test_that("out-of-sample prediction returns nonnegative weights per level (K=3)", {
  d <- make_data()
  fit <- energy_balance(x = d$X, treatment = d$tm)

  set.seed(20260811L)
  m <- 40
  Xnew <- matrix(rnorm(m * d$p), m, d$p)
  Anew <- sample(1:3, m, replace = TRUE)
  wp <- predict(fit, newdata = Xnew, treatment = Anew)

  expect_length(wp, m)
  expect_true(all(wp >= 0))
  for (lv in 1:3)
    expect_true(all(wp[Anew == lv] >= 0))

  # a treatment level not seen at fit time is an error
  expect_error(predict(fit, newdata = Xnew, treatment = rep(9L, m)),
               "among the fitted levels")
})

test_that("the formula interface routes a 3-level factor to the grouped engine", {
  set.seed(20260810L)
  n <- 240; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  df <- as.data.frame(X); names(df) <- paste0("V", seq_len(p))
  df$A <- factor(sample(c("a", "b", "c"), n, replace = TRUE))

  fit <- energy_balance(A ~ ., data = df)
  expect_s3_class(fit, "ebw")
  expect_equal(fit$treatment_type, "multinomial")
  expect_equal(fit$n_groups, 3L)

  # formula-carried treatment lets predict map groups out of sample
  set.seed(20260811L)
  nd <- as.data.frame(matrix(rnorm(20 * p), 20, p))
  names(nd) <- paste0("V", seq_len(p))
  nd$A <- factor(sample(c("a", "b", "c"), 20, replace = TRUE), levels = c("a", "b", "c"))
  wp <- predict(fit, newdata = nd)
  expect_length(wp, 20)
  expect_true(all(wp >= 0))
})
