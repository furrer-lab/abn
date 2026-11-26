# Scenario: a directed graph A -> B -> C with an extra A -> C edge.
# Inputs: labeled igraph with three nodes.
# Expected: parents(C) = {A, B}; children(A) = {B, C}; roots and leaves are empty.
test_that("find_parents() and find_children() return exact neighbors", {
  graph <- igraph::graph_from_literal(A -+ B, B -+ C, A -+ C)

  expect_setequal(find_parents(graph, "C"), c("A", "B"))
  expect_setequal(find_children(graph, "A"), c("B", "C"))
  expect_length(find_parents(graph, "A"), 0)
  expect_length(find_children(graph, "C"), 0)
})

# Scenario: unlabeled, out-of-range, unknown, and undirected graph inputs.
# Inputs: igraph objects without names, a missing vertex label, and an undirected graph.
# Expected: documented errors for invalid inputs and a warning plus neighbor output for the undirected graph.
test_that("find_parents() and find_children() handle invalid graph inputs", {
  unlabeled <- igraph::make_graph(c(1, 2), directed = TRUE)
  expect_error(find_parents(unlabeled, "A"), "not labeled")
  expect_error(find_parents(unlabeled, 3), "less than")

  labeled <- igraph::graph_from_literal(A -+ B)
  expect_error(find_parents(labeled, "Z"), "not part of the graph")

  undirected <- igraph::make_graph(c(1, 2), directed = FALSE)
  expect_warning(parents <- find_parents(undirected, 1), "not directed")
  expect_setequal(parents, "2")
})

# Scenario: valid binomial evidence for a node with levels n/y.
# Inputs: a factor column and an evidence value that belongs to its levels.
# Expected: the returned evidence is a factor with the data's levels and the supplied value.
test_that("check_evidence() preserves valid binomial evidence as a factor", {
  data <- data.frame(
    b1 = factor(c("n", "y"), levels = c("n", "y")),
    p1 = c(1, 2),
    g1 = c(3, 4)
  )
  dists <- list(b1 = "binomial", p1 = "poisson", g1 = "gaussian")

  evidence <- check_evidence(data, dists, hypothesis = NULL,
                             evidence = list(b1 = "y"))

  expect_s3_class(evidence$b1, "factor")
  expect_equal(levels(evidence$b1), c("n", "y"))
  expect_equal(as.character(evidence$b1), "y")
})

# Scenario: invalid evidence values, an unknown node, and a character value for a numeric node.
# Inputs: evidence entries that violate the distribution contract.
# Expected: invalid entries are discarded, with a specific warning plus a summary
#           "wrong information" warning.
test_that("check_evidence() discards invalid evidence with warnings", {
  data <- data.frame(
    b1 = factor(c("n", "y"), levels = c("n", "y")),
    p1 = c(1, 2),
    g1 = c(3, 4)
  )
  dists <- list(b1 = "binomial", p1 = "poisson", g1 = "gaussian")

  collect_warnings <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(value = value, warnings = warnings)
  }

  res <- collect_warnings(
    check_evidence(data, dists, NULL, list(b1 = "z"))
  )
  expect_length(res$value, 0)
  expect_true(any(grepl("does not have an expected value", res$warnings)))
  expect_true(any(grepl("wrong information", res$warnings)))

  res <- collect_warnings(
    check_evidence(data, dists, NULL, list(nope = "x"))
  )
  expect_length(res$value, 0)
  expect_true(any(grepl("wrong information", res$warnings)))

  res <- collect_warnings(
    check_evidence(data, dists, NULL, list(p1 = "3"))
  )
  expect_length(res$value, 0)
  expect_true(any(grepl("should be numeric", res$warnings)))
  expect_true(any(grepl("wrong information", res$warnings)))
})

# Scenario: the hypothesis node is also supplied as evidence.
# Inputs: hypothesis b1 and evidence b1 = "y".
# Expected: a printed notice and only b1 returned, unchanged.
test_that("check_evidence() keeps only the hypothesis when it is evidence", {
  data <- data.frame(b1 = factor(c("n", "y"), levels = c("n", "y")))
  dists <- list(b1 = "binomial")

  expect_output(
    ev <- check_evidence(data, dists, hypothesis = "b1",
                         evidence = list(b1 = "y")),
    "already known"
  )

  expect_equal(ev, list(b1 = "y"))
})

# Scenario: Gaussian root with a known intercept and variance.
# Inputs: intercept 2, precision coefficient 0.25, variance 4; also var = NULL.
# Expected: c(mean, variance) = c(2, 4), and c(2, 1) when var is absent.
test_that("predict_root() returns the Gaussian mean and variance", {
  fit <- list(g1 = c("g1|intercept" = 2, "precision" = 0.25))
  data <- data.frame(g1 = 0)
  dists <- list(g1 = "gaussian")

  expect_equal(
    unname(predict_root(fit, c(g1 = 4), data, dists, "g1")),
    c(2, 4)
  )
  expect_equal(
    unname(predict_root(fit, NULL, data, dists, "g1")),
    c(2, 1)
  )
})

# Scenario: Poisson root with a known log intercept.
# Inputs: intercept log(3).
# Expected: rate exp(log(3)) = 3.
test_that("predict_root() returns the Poisson rate", {
  fit <- list(p1 = c("p1|intercept" = log(3)))
  data <- data.frame(p1 = 0)
  dists <- list(p1 = "poisson")

  expect_equal(unname(predict_root(fit, NULL, data, dists, "p1")), 3)
})

# Scenario: binomial root with a known logit intercept.
# Inputs: intercept qlogis(0.7) and levels n/y.
# Expected: probabilities n = 0.3 and y = 0.7.
test_that("predict_root() returns the binomial probabilities", {
  fit <- list(b1 = c("b1|intercept" = stats::qlogis(0.7)))
  data <- data.frame(b1 = factor(c("n", "y"), levels = c("n", "y")))
  dists <- list(b1 = "binomial")

  expect_equal(
    predict_root(fit, NULL, data, dists, "b1"),
    c(n = 1 - 0.7, y = 0.7)
  )
})

# Scenario: multinomial root with known category logits and a reference category.
# Inputs: logits log(2) and log(3) for categories 2 and 3; category 1 is the reference.
# Expected: softmax probabilities 1/6, 2/6, 3/6 for categories 1, 2, 3.
test_that("predict_root() returns the multinomial softmax probabilities", {
  fit <- list(
    m1 = c("m1|intercept.2" = log(2), "m1|intercept.3" = log(3))
  )
  data <- data.frame(m1 = factor(c("1", "2", "3"), levels = c("1", "2", "3")))
  dists <- list(m1 = "multinomial")

  expect_equal(
    predict_root(fit, NULL, data, dists, "m1"),
    c("1" = 1 / 6, "2" = 2 / 6, "3" = 3 / 6)
  )
})

# Scenario: a root distribution not supported by predict_root().
# Inputs: distribution "weibull".
# Expected: an error naming the unsupported distribution.
test_that("predict_root() rejects unsupported distributions", {
  expect_error(
    predict_root(
      list(x = c("x|intercept" = 1)),
      NULL,
      data.frame(x = 0),
      list(x = "weibull"),
      "x"
    ),
    "Unsupported distribution type"
  )
})

# Scenario: Poisson root and one-parent Poisson node with known coefficients.
# Inputs: root intercept log(3); child intercept log(2), slope 0.5, parent value 2.
# Expected: root rate 3; child rate exp(log(2) + 0.5 * 2).
test_that("predict_node_from_parent_poisson() uses the exact log link", {
  data <- data.frame(p1 = 2, p2 = 0)
  dists <- list(p1 = "poisson", p2 = "poisson")

  root_fit <- list(p1 = c("p1|intercept" = log(3)))
  expect_equal(
    unname(predict_node_from_parent_poisson(
      data, dists, root_fit, "p1", list(), character(0)
    )),
    3
  )

  child_fit <- list(p2 = c("p2|intercept" = log(2), p1 = 0.5))
  expect_equal(
    unname(predict_node_from_parent_poisson(
      data, dists, child_fit, "p2", list(), "p1", list(p1 = 2)
    )),
    exp(log(2) + 0.5 * 2)
  )
})

# Scenario: Gaussian child with one parent and a known variance, plus evidence passthrough.
# Inputs: intercept 1, slope 2, parent value 3, variance 9; evidence g2 = 5.
# Expected: c(7, 9) for the model; c(5, 1) when g2 is evidence.
test_that("predict_node_from_parent_gaussian() returns known mean and variance", {
  data <- data.frame(g1 = 3, g2 = 0)
  dists <- list(g1 = "gaussian", g2 = "gaussian")
  fit <- list(g2 = c("g2|intercept" = 1, g1 = 2))
  var <- c(g2 = 9)

  expect_equal(
    unname(predict_node_from_parent_gaussian(
      data, dists, fit, var, "g2", list(), "g1", list(g1 = 3)
    )),
    c(7, 9)
  )

  expect_equal(
    unname(predict_node_from_parent_gaussian(
      data, dists, fit, var, "g2", list(g2 = 5), "g1", list(g1 = 3)
    )),
    c(5, 1)
  )
})

# Scenario: binomial child with one Gaussian parent and known coefficients.
# Inputs: intercept qlogis(0.2), slope 1, parent value 0, levels n/y.
# Expected: probabilities n = 0.8 and y = 0.2.
test_that("predict_node_from_parent_binomial() returns known probabilities", {
  data <- data.frame(
    g1 = 0,
    b2 = factor(c("n", "y"), levels = c("n", "y"))
  )
  dists <- list(g1 = "gaussian", b2 = "binomial")
  fit <- list(b2 = c("b2|intercept" = stats::qlogis(0.2), g1 = 1))

  expect_equal(
    predict_node_from_parent_binomial(
      data, dists, fit, "b2", list(), "g1", list(g1 = 0)
    ),
    c(n = 1 - 0.2, y = 0.2)
  )
})

# Scenario: multinomial child with one Gaussian parent and known category coefficients.
# Inputs: child levels 1/2/3; intercepts log(2), log(3); slopes 0.5, -0.5; parent value 2.
# Expected: normalized probabilities 1/D, 2*exp(1)/D, 3*exp(-1)/D
#           with D = 1 + 2*exp(1) + 3*exp(-1).
test_that("predict_node_from_parent_multinomial() returns known softmax probabilities", {
  data <- data.frame(
    m2 = factor(c("1", "2", "3"), levels = c("1", "2", "3")),
    g1 = 2
  )
  dists <- list(m2 = "multinomial", g1 = "gaussian")
  fit <- list(
    m2 = c(
      "m2|intercept.2" = log(2),
      "m2|intercept.3" = log(3),
      "g1.2" = 0.5,
      "g1.3" = -0.5
    )
  )

  expected <- c(
    "1" = 1,
    "2" = 2 * exp(1),
    "3" = 3 * exp(-1)
  )
  expected <- expected / sum(expected)

  expect_equal(
    predict_node_from_parent_multinomial(
      data, dists, fit, "m2", list(), "g1", list(g1 = 2)
    ),
    expected
  )
})

# Scenario: Poisson prior and log-likelihood against base R.
# Inputs: counts 0, 1, 2; rate 3; y = 2, x = 1.5, coef = 0.5, continuous_part = 0.2.
# Expected: exact dpois values and the dpois log-density.
test_that("Poisson prior and log-likelihood match base R", {
  expect_equal(
    prior_poisson(c(0, 1, 2), 3),
    stats::dpois(c(0, 1, 2), 3)
  )
  expect_equal(
    LogL_poisson(2, 1.5, 0.5, 0.2),
    stats::dpois(2, exp(0.5 * 1.5 + 0.2), log = TRUE),
    tolerance = 1e-12
  )
})

# Scenario: Gaussian prior and log-likelihood against base R.
# Inputs: x = 0, mu = 1, variance 4; and y = 1, x = 2, coef = 0.5, var = 4, continuous_part = -1.
# Expected: exact dnorm density and dnorm log-density.
test_that("Gaussian prior and log-likelihood match base R", {
  expect_equal(
    prior_gaussian(0, 1, 4),
    stats::dnorm(0, mean = 1, sd = 2)
  )
  expect_equal(
    LogL_gaussian(1, 2, 0.5, 4, -1),
    stats::dnorm(1, mean = -1 + 0.5 * 2, sd = 2, log = TRUE),
    tolerance = 1e-12
  )
})

# Scenario: binomial prior and log-likelihood against base R, including probability clamping.
# Inputs: x = 1 with p = 0.7, p = 1.2, p = -0.2; and y = 1, x = 2, coef = 0.5, continuous_part = -1.
# Expected: exact dbinom probabilities after clamping and the dbinom log-density.
test_that("Binomial prior and log-likelihood match base R", {
  expect_equal(prior_binomial(1, 0.7), stats::dbinom(1, 1, 0.7))
  expect_equal(prior_binomial(1, 1.2), stats::dbinom(1, 1, 1))
  expect_equal(prior_binomial(0, -0.2), stats::dbinom(0, 1, 0))
  expect_equal(
    LogL_binomial(1, 2, 0.5, -1),
    stats::dbinom(1, 1, stats::plogis(-1 + 0.5 * 2), log = TRUE),
    tolerance = 1e-12
  )
})

# Scenario: multinomial log-likelihood corresponds to the softmax probability of the observed category.
# Inputs: three categories, coefficients (0.5, -0.2), intercepts (0.1, 0.3), parent value 2.
# Expected: exp(logL_multinomial(y)) equals the softmax probability for y, for y = 1, 2, 3.
test_that("logL_multinomial() matches the softmax probabilities", {
  coef <- c(0.5, -0.2)
  continuous_part <- c(0.1, 0.3)
  x <- 2

  eta <- c(0, continuous_part[1] + coef[1] * x,
           continuous_part[2] + coef[2] * x)
  probs <- exp(eta) / sum(exp(eta))

  for (y in 1:3) {
    expect_equal(
      exp(logL_multinomial(y, x, coef, continuous_part)),
      probs[y]
    )
  }
})
