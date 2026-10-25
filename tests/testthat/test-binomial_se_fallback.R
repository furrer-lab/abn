test_that("SVD fallback matches solve() on well-conditioned PSD matrices", {
  se_normal_path <- function(M) {
    as.matrix(sqrt(diag(solve(M))))
  }
  se_svd_fallback <- function(M) {
    s <- svd(M)
    tol <- max(s$d) * .Machine$double.eps * length(s$d)
    d_inv <- ifelse(s$d > tol, 1 / s$d, 0)
    as.matrix(sqrt(rowSums(s$v^2 * rep(d_inv, each = nrow(s$v)))))
  }
  make_psd <- function(eigvals, seed = 1) {
    set.seed(seed)
    p <- length(eigvals)
    V <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
    V %*% diag(eigvals) %*% t(V)
  }

  cases <- list(
    identity = make_psd(rep(1, 5), seed = 2),
    moderate = make_psd(c(1000, 100, 10, 5, 1), seed = 3),
    larger_rng = make_psd(sort(runif(12, 0.5, 20), decreasing = TRUE), seed = 4)
  )
  for (nm in names(cases)) {
    M <- unname(cases[[nm]])
    expect_equal(se_svd_fallback(M), se_normal_path(M), tolerance = 1e-8, info = nm)
  }
})

test_that("SVD fallback matches solve() on a realistic logistic XtWX", {
  se_normal_path <- function(M) {
    as.matrix(sqrt(diag(solve(M))))
  }
  se_svd_fallback <- function(M) {
    s <- svd(M)
    tol <- max(s$d) * .Machine$double.eps * length(s$d)
    d_inv <- ifelse(s$d > tol, 1 / s$d, 0)
    as.matrix(sqrt(rowSums(s$v^2 * rep(d_inv, each = nrow(s$v)))))
  }

  set.seed(7)
  n <- 200
  X <- cbind(1, matrix(rnorm(n * 3), n, 3))
  mu <- plogis(X %*% c(-0.3, 0.5, -0.2, 0.4))
  M <- crossprod(X, as.numeric(mu * (1 - mu)) * X)
  expect_equal(se_svd_fallback(M), se_normal_path(M), tolerance = 1e-8)
})

test_that("SVD fallback matches solve() on a varcov shaped like a real abn binomial node", {
  se_normal_path <- function(M) {
    as.matrix(sqrt(diag(solve(M))))
  }
  se_svd_fallback <- function(M) {
    s <- svd(M)
    tol <- max(s$d) * .Machine$double.eps * length(s$d)
    d_inv <- ifelse(s$d > tol, 1 / s$d, 0)
    as.matrix(sqrt(rowSums(s$v^2 * rep(d_inv, each = nrow(s$v)))))
  }

  mydat <- ex0.dag.data[, c("b1", "b2", "g1")]
  glm_fit <- glm(b1 ~ b2 + g1, data = mydat, family = binomial())
  XtWX <- unname(solve(vcov(glm_fit)))
  expect_equal(se_svd_fallback(XtWX), se_normal_path(XtWX), tolerance = 1e-8)
})

test_that("fitAbn binomial fallback returns finite SEs on near-singular varcov", {
  skip_on_cran()
  set.seed(1)
  n <- 50
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(2 * x))
  x2 <- x + rnorm(n, sd = 1e-10)
  df <- data.frame(y = as.factor(y), x = x, x2 = x2)
  dists <- list(y = "binomial", x = "gaussian", x2 = "gaussian")
  dag <- matrix(0, 3, 3, dimnames = list(names(df), names(df)))
  dag["y", c("x", "x2")] <- 1
  suppressMessages({
    fit <- fitAbn(dag = dag, data.df = df, data.dists = dists, method = "mle")
  })
  se <- fit$Stderror[["y"]]
  se_nonna <- se[!is.na(se)]
  expect_true(length(se_nonna) > 0)
  expect_true(all(is.finite(se_nonna)))
  expect_true(all(se_nonna > 0))
})
