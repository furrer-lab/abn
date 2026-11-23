# Synthetic grouped data for tests that need all four distributions.
#
# The fixture is generated here, deterministically, instead of being stored as
# a binary .rda. It has a real group identifier with 20 groups and 20
# observations per group. In particular, this avoids using FCV$Pedigree as a
# grouping variable: Pedigree is a binary cat-level trait, not a cluster ID,
# and the resulting two-group random-effects multinomial fits are unstable.
#
# Generative model:
#   u_g ~ N(0, 0.5^2)                            group-level random intercept
#   x   ~ N(0, 1)                                latent covariate
#   g   = 0.5 + 0.4 x + u_g + N(0, 0.5^2)        Gaussian node
#   p   ~ Poisson(exp(0.5 + 0.2 x + u_g))        count node
#   b   ~ Bernoulli(plogis(-0.2 + 0.4 x + u_g))  binary node
#   m   ~ Categorical(softmax(0.4 + 0.3 x + u_g,
#                             -0.3 + 0.2 x + u_g, 0))
#                                                 multinomial node
#
# Test DAG: p <- g, b <- g, m is a root.
#
# Keeping m as a root is deliberate. It exercises grouped multinomial fitting
# and simulation without relying on a near-singular parent model. The seed is
# fixed so the fixture is reproducible.
make_grouped_test_data <- function(seed = 20261123,
                                   n_groups = 20,
                                   n_per_group = 20) {
  # Do not disturb the RNG state of the calling test.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  n <- n_groups * n_per_group
  group <- factor(rep(seq_len(n_groups), each = n_per_group))
  u <- rnorm(n_groups, sd = 0.5)[group]
  x <- rnorm(n)

  g <- 0.5 + 0.4 * x + u + rnorm(n, sd = 0.5)
  p <- rpois(n, exp(0.5 + 0.2 * x + u))
  b <- factor(rbinom(n, 1, plogis(-0.2 + 0.4 * x + u)),
              levels = c(0, 1))

  eta1 <- 0.4 + 0.3 * x + u
  eta2 <- -0.3 + 0.2 * x + u
  probs <- cbind(exp(eta1), exp(eta2), 1)
  probs <- probs / rowSums(probs)
  m <- factor(
    apply(probs, 1, function(pr) sample.int(3, 1, prob = pr)),
    levels = 1:3
  )

  stopifnot(
    nlevels(droplevels(b)) == 2,
    nlevels(droplevels(m)) == 3
  )

  data <- data.frame(group = group, g = g, p = p, b = b, m = m)
  dists <- list(g = "gaussian", p = "poisson",
                b = "binomial", m = "multinomial")
  dag <- matrix(0, 4, 4, dimnames = list(names(dists), names(dists)))
  dag["p", "g"] <- 1
  dag["b", "g"] <- 1

  list(data = data, dists = dists, dag = dag, group.var = "group")
}
