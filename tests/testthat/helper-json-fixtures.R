json_fixture_empty_dag <- function(dists) {
  dag <- matrix(0, nrow = length(dists), ncol = length(dists))
  colnames(dag) <- rownames(dag) <- names(dists)
  dag
}

json_fixture_assert_abnfit <- function(fit, method, dists, grouped = FALSE,
                                       multinomial = FALSE) {
  expect_s3_class(fit, "abnFit")
  expect_s3_class(fit$abnDag, "abnDag")
  expect_equal(fit$method, method)
  expect_named(fit$abnDag$data.dists, names(dists))
  expect_equal(unname(unlist(fit$abnDag$data.dists)), unname(unlist(dists)))
  expect_equal(dim(fit$abnDag$dag), c(length(dists), length(dists)))
  expect_setequal(colnames(fit$abnDag$dag), names(dists))
  expect_setequal(rownames(fit$abnDag$dag), names(dists))

  if (!grouped) {
    expect_equal(length(fit$coef), length(dists))
    expect_named(fit$coef, names(dists))
    expect_named(fit$Stderror, names(dists))
  } else {
    if (!is.null(fit$coef)) expect_named(fit$coef, names(dists))
    if (!is.null(fit$Stderror)) expect_named(fit$Stderror, names(dists))
  }

  if (multinomial) {
    expect_true("multinomial" %in% unlist(dists))
    expect_true("multinomial" %in% unlist(fit$abnDag$data.dists))
    expect_true(!is.null(fit$multinomial.states) ||
                  any(unlist(fit$abnDag$data.dists) == "multinomial"))
  }

  if (grouped) {
    expect_true(!is.null(fit$group.var) || !is.null(fit$mu))
    expect_true(!is.null(fit$mu) || !is.null(fit$betas) ||
                  !is.null(fit$sigma) || !is.null(fit$sigma_alpha))
  }

  invisible(fit)
}

json_fixture_fit_g2b2c_mle <- function() {
  dists <- list(
    G1 = "gaussian",
    B1 = "binomial",
    B2 = "binomial",
    C = "multinomial",
    G2 = "gaussian"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["B1", "G1"] <- 1
  dag["B2", "B1"] <- 1
  dag["C", c("B1", "B2")] <- 1
  dag["G2", c("G1", "C")] <- 1

  fitAbn(dag = dag, data.df = g2b2c_data, data.dists = dists,
         method = "mle", centre = FALSE)
}

json_fixture_fit_g2pbcgrp_mle_grouped <- function() {
  dists <- list(
    G1 = "gaussian",
    P = "poisson",
    B = "binomial",
    C = "multinomial",
    G2 = "gaussian"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["P", "G1"] <- 1
  dag["B", "G1"] <- 1
  dag["C", c("P", "B")] <- 1
  dag["G2", c("G1", "C")] <- 1

  fitAbn(dag = dag, data.df = droplevels(g2pbcgrp), data.dists = dists,
         method = "mle", group.var = "group", centre = FALSE)
}

json_fixture_fit_ex1_mle <- function() {
  data.df <- droplevels(ex1.dag.data[seq_len(1000),
                                    c("b1", "p1", "g1", "b2", "p2", "b3", "g2")])
  dists <- list(
    b1 = "binomial",
    p1 = "poisson",
    g1 = "gaussian",
    b2 = "binomial",
    p2 = "poisson",
    b3 = "binomial",
    g2 = "gaussian"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["p2", c("b1", "p1")] <- 1
  dag["b3", c("b1", "g1", "b2")] <- 1
  dag["g2", c("p1", "g1", "b2")] <- 1

  fitAbn(dag = dag, data.df = data.df, data.dists = dists,
         method = "mle", centre = FALSE)
}

json_fixture_fit_ex1_bayes <- function() {
  data.df <- droplevels(ex1.dag.data[seq_len(1000), c("b1", "p1", "g1", "p2")])
  dists <- list(
    b1 = "binomial",
    p1 = "poisson",
    g1 = "gaussian",
    p2 = "poisson"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["p2", c("b1", "p1")] <- 1

  fitAbn(dag = dag, data.df = data.df, data.dists = dists,
         method = "bayes", centre = FALSE)
}

json_fixture_fit_ex3_mle_grouped <- function() {
  dists <- list(b1 = "binomial", b2 = "binomial")

  fitAbn(dag = ~ b1 | b2, data.df = ex3.dag.data[, c(1, 2, 14)],
         data.dists = dists, method = "mle", group.var = "group",
         cor.vars = c("b1", "b2"), centre = FALSE)
}

json_fixture_fit_adg_mle_grouped <- function() {
  data.df <- droplevels(adg[, c("AR", "eggs", "wormCount", "age", "adg", "farm")])
  data.df$AR <- factor(data.df$AR)
  data.df$eggs <- factor(data.df$eggs)
  data.df$farm <- factor(data.df$farm)

  dists <- list(
    AR = "binomial",
    eggs = "binomial",
    wormCount = "poisson",
    age = "gaussian",
    adg = "gaussian"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["eggs", "AR"] <- 1
  dag["wormCount", "eggs"] <- 1
  dag["adg", c("age", "wormCount")] <- 1

  fitAbn(dag = dag, data.df = data.df, data.dists = dists,
         method = "mle", group.var = "farm", centre = FALSE)
}

json_fixture_fit_fcv_mle <- function() {
  data.df <- droplevels(FCV[, c("Outdoor", "Sex", "GroupSize", "Age")])
  dists <- list(
    Outdoor = "binomial",
    Sex = "multinomial",
    GroupSize = "poisson",
    Age = "gaussian"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["Outdoor", "Sex"] <- 1
  dag["GroupSize", "Outdoor"] <- 1
  dag["Age", c("Sex", "GroupSize")] <- 1

  fitAbn(dag = dag, data.df = data.df, data.dists = dists,
         method = "mle", centre = FALSE)
}
