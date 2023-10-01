test_that("Gaussian response Bugs model construction works", {
  if(.Platform$OS.type == "unix") {
    capture.output(
      # No parent nodes / no predictors
      expect_no_error(gauss_bugs(nodename = "a",
                                 parentnames = NULL,
                                 nodesintercept = c(0.318077),
                                 parentcoefs = NULL,
                                 std = c(0.05773503))),

      # One parent nodes / predictor
      expect_no_error(gauss_bugs(nodename = "a",
                                 parentnames = "b",
                                 nodesintercept = c(0.318077),
                                 parentcoefs = list("b"=c(b=0.3059395)),
                                 std = c(0.05773503))),

      # Multiple parent nodes / predictors
      expect_no_error(gauss_bugs(nodename = "a",
                                 parentnames = c("b", "c"),
                                 nodesintercept = c(0.318077),
                                 parentcoefs = list("b"=c(b=0.3059395),
                                                    "c"=c(c=0.5555)),
                                 std = c(0.05773503))),
      # Multinomial parent nodes / predictors
      expect_no_error(gauss_bugs(nodename = "a",
                                 parentnames = c("b", "c"),
                                 nodesintercept = c(0.318077),
                                 parentcoefs = list("b"=c(b=0.3059395, b=0.66666, b=0.77777), # multinomial with three factor levels
                                                    "c"=c(c=0.5555)),
                                 std = c(0.05773503))),
      file = "/dev/null"
    )
  }

  outstring <- c(
    "a ~ dnorm(mu.a, precision.a) # Gaussian response
mu.a <- 0.318077 + 0.3059395*b + 0.5555*c # Linear regression
precision.a <- inverse(0.05773503) # precision tau = 1/standard_dev")
  expect_output(gauss_bugs(nodename = "a",
                           parentnames = c("b", "c"),
                           nodesintercept = c(0.318077),
                           parentcoefs = list("b"=c(b=0.3059395),
                                              "c"=c(c=0.5555)),
                           std = c(0.05773503)),
                regexp = outstring,
                fixed = TRUE)
})

test_that("Gaussian response Bugs model construction with mixed-effects works", {
  muY <- 3
  sigmaY <- 1
  sigma_alphaY <- 0.5
  betaX <- 4
  betaN <- 5
  if(.Platform$OS.type == "unix") {
    capture.output(
      # No parent nodes / no predictors
      expect_no_error(gauss_bugsGroup(nodename = "Y",
                                      nodesintercept = muY,
                                      parentnames = NULL,
                                      parentcoefs = NULL,
                                      sigma = sigmaY,
                                      sigma_alpha = sigma_alphaY)),

      # One parent nodes / predictor
      expect_no_error(gauss_bugsGroup(nodename = "Y",
                                      nodesintercept = muY,
                                      parentnames = c("X"),
                                      parentcoefs = list("X"=c(X=betaX)),
                                      sigma = sigmaY,
                                      sigma_alpha = sigma_alphaY)),

      # Multiple parent nodes / predictors
      expect_no_error(gauss_bugsGroup(nodename = "Y",
                                      nodesintercept = muY,
                                      parentnames = c("X", "N"),
                                      parentcoefs = list("X"=c(X=betaX),
                                                         "N"=c(N=betaN)),
                                      sigma = sigmaY,
                                      sigma_alpha = sigma_alphaY)),

      # Multinomial parent nodes / predictors
      # expect_no_error(),
      file = "/dev/null"
    )
  }

  outstring <- c(
    "Y <- mu_Y + 4*X + 5*N + alpha_Y + e_Y
mu_Y <- 3
e_Y ~ dnorm(mu_e_Y, tau_Y)
mu_e_Y <- 0
tau_Y <- inverse(sigma_Y)
sigma_Y <- 1
alpha_Y ~ dnorm(mu_alpha_Y, tau_alpha_Y)
mu_alpha_Y <- 0
tau_alpha_Y <- inverse(sigma_alpha_Y)
sigma_alpha_Y <- 0.5")
  expect_output(gauss_bugsGroup(nodename = "Y",
                                nodesintercept = muY,
                                parentnames = c("X", "N"),
                                parentcoefs = list("X"=c(X=betaX),
                                                   "N"=c(N=betaN)),
                                sigma = sigmaY,
                                sigma_alpha = sigma_alphaY),
                regexp = outstring,
                fixed = TRUE)
})

test_that("Binomial response Bugs model construction works", {
  if(.Platform$OS.type == "unix") {
    capture.output(
      # No parent nodes / no predictors
      expect_no_error(bern_bugs(nodename = "a",
                                parentnames = NULL,
                                nodesintercept = c(0.318077),
                                parentcoefs = NULL)),

      # One parent nodes / predictor
      expect_no_error(bern_bugs(nodename = "a",
                                parentnames = "b",
                                nodesintercept = c(0.318077),
                                parentcoefs = list("b"=c(b=0.3059395)))),

      # Multiple parent nodes / predictors
      expect_no_error(bern_bugs(nodename = "a",
                                parentnames = c("b", "c"),
                                nodesintercept = c(0.318077),
                                parentcoefs = list("b"=c(b=0.3059395),
                                                   "c"=c(c=0.5555)))),
      file = "/dev/null"
    )
  }

  outstring <- c(
    "a ~ dbern(p.a) # Bernoulli response
logit(p.a) <- 0.318077 + 0.3059395*b + 0.5555*c # logistic regression")
  expect_output(bern_bugs(nodename = "a",
                          parentnames = c("b", "c"),
                          nodesintercept = c(0.318077),
                          parentcoefs = list("b"=c(b=0.3059395),
                                             "c"=c(c=0.5555))),
                regexp = outstring,
                fixed = TRUE)
})

test_that("Binomial response Bugs model construction with mixed-effects works", {
  muY <- 3
  sigma_alphaY <- 0.5
  betaX <- 4
  betaN <- 5
  if(.Platform$OS.type == "unix") {
    capture.output(
      # No parent nodes / no predictors
      expect_no_error(bern_bugsGroup(nodename = "Y",
                                     nodesintercept = muY,
                                     parentnames = NULL,
                                     parentcoefs = NULL,
                                     sigma_alpha = sigma_alphaY)),

      # One parent nodes / predictor
      expect_no_error(bern_bugsGroup(nodename = "Y",
                                     nodesintercept = muY,
                                     parentnames = c("X"),
                                     parentcoefs = list("X"=c(X=betaX)),
                                     sigma_alpha = sigma_alphaY)),

      # Multiple parent nodes / predictors
      expect_no_error(bern_bugsGroup(nodename = "Y",
                                     nodesintercept = muY,
                                     parentnames = c("X", "N"),
                                     parentcoefs = list("X"=c(X=betaX),
                                                        "N"=c(N=betaN)),
                                     sigma_alpha = sigma_alphaY)),

      # Multinomial parent nodes / predictors
      # expect_no_error(),
      file = "/dev/null"
    )
  }

  outstring <- c(
    "Y ~ dbern(p_Y)
logit(p_Y) <- mu_Y + 4*X + 5*N + alpha_Y
mu_Y <- 3
alpha_Y ~ dnorm(mu_alpha_Y, tau_alpha_Y)
mu_alpha_Y <- 0
tau_alpha_Y <- inverse(sigma_alpha_Y)
sigma_alpha_Y <- 0.5")
  expect_output(bern_bugsGroup(nodename = "Y",
                               nodesintercept = muY,
                               parentnames = c("X", "N"),
                               parentcoefs = list("X"=c(X=betaX),
                                                  "N"=c(N=betaN)),
                               sigma_alpha = sigma_alphaY),
                regexp = outstring,
                fixed = TRUE)
})

test_that("Categorical response Bugs model construction works", {
  if(.Platform$OS.type == "unix") {
    capture.output(
      # No parent nodes / no predictors
      expect_no_error(categorical_bugs(nodename = "b",
                                       nodesCatIdx = c(2, 3, 4),
                                       parentnames = NULL,
                                       nodesintercepts = c(2.188650, 3.133928, 3.138531),
                                       parentcoefs = NULL)),

      # One parent nodes / predictor
      expect_no_error(categorical_bugs(nodename = "b",
                                       nodesCatIdx = c(2, 3, 4),
                                       parentnames = "a",
                                       nodesintercepts = c(2.188650, 3.133928, 3.138531),
                                       parentcoefs = list("a"=c(a=1.686432, a=3.134161, a=5.052104)))),

      # Multiple parent nodes / predictors
      expect_no_error(categorical_bugs(nodename = "b",
                                       nodesCatIdx = c(2, 3, 4),
                                       parentnames = c("a", "c"),
                                       nodesintercepts = c(2.188650, 3.133928, 3.138531),
                                       parentcoefs = list("a"=c(a=1.686432, a=3.134161, a=5.052104),
                                                          "c"=c(c=0.5555, c=0.6666, c=0.7777)))),
      file = "/dev/null"
    )
  }

  outstring <- c(
    "b ~ dcat(p.b) # Categorical response
p.b[1] <- phi.b[1]/sum(phi.b) # soft-max
log(phi.b[1]) <- 0 # Reference category
p.b[2] <- phi.b[2]/sum(phi.b) # soft-max
log(phi.b[2]) <- 2.18865 + 1.686432*a
p.b[3] <- phi.b[3]/sum(phi.b) # soft-max
log(phi.b[3]) <- 3.133928 + 3.134161*a
p.b[4] <- phi.b[4]/sum(phi.b) # soft-max
log(phi.b[4]) <- 3.138531 + 5.052104*a")
  expect_output(categorical_bugs(nodename = "b",
                                 nodesCatIdx = c(2, 3, 4),
                                 parentnames = "a",
                                 nodesintercepts = c(2.188650, 3.133928, 3.138531),
                                 parentcoefs = list("a"=c(a=1.686432, a=3.134161, a=5.052104))),
                regexp = outstring,
                fixed = TRUE)
})

test_that("Categorical response Bugs model construction with mixed-effects works", {
  nodename <- "Y"
  nodesCatIdx <- c(2,3,4)
  nodesintercepts <- c(5,6,7)
  parentnames <- c("X")
  parentcoefs <- data.frame(X = c(4,5,6))
  sigma <- 1 # within variance
  sigma_alpha <- matrix(0, ncol = 4, nrow = 4) # var-covariance for one predictor
  diag(sigma_alpha) <- 0.5
  sigma_alpha2 <- matrix(0, ncol = 8, nrow = 8) # var-covariance for two predictors
  diag(sigma_alpha2) <- 0.5
  if(.Platform$OS.type == "unix") {
    capture.output(
      # No parent nodes / no predictors
      expect_no_error(categorical_bugsGroup(nodename = nodename,
                                            nodesCatIdx = nodesCatIdx,
                                            nodesintercepts = nodesintercepts,
                                            parentnames = NULL,
                                            parentcoefs = NULL,
                                            sigma = sigma,
                                            sigma_alpha = sigma_alpha)),

      # One parent nodes / predictor
      expect_no_error(categorical_bugsGroup(nodename = nodename,
                                            nodesCatIdx = nodesCatIdx,
                                            nodesintercepts = nodesintercepts,
                                            parentnames = parentnames,
                                            parentcoefs = parentcoefs,
                                            sigma = sigma,
                                            sigma_alpha = sigma_alpha)),

      # Multiple parent nodes / predictors
      expect_no_error(categorical_bugsGroup(nodename = nodename,
                                            nodesCatIdx = nodesCatIdx,
                                            nodesintercepts = nodesintercepts,
                                            parentnames = c("X1", "X2"),
                                            parentcoefs =  data.frame(X1 = c(4,5,6),
                                                                      X2 = c(4,5,6)),
                                            sigma = sigma,
                                            sigma_alpha = sigma_alpha2)),

      # Multinomial parent nodes / predictors
      # expect_no_error(),
      file = "/dev/null"
    )
  }

  outstring <- c(
    "Y ~ dcat(p_Y)
p_Y[1] <- phi_Y[1]/sum(phi_Y)
log(phi_Y[1]) <- 0 + alpha_Y[1]
p_Y[2] <- phi_Y[2]/sum(phi_Y)
log(phi_Y[2]) <- 5 + 4*X + alpha_Y[1]
p_Y[3] <- phi_Y[3]/sum(phi_Y)
log(phi_Y[3]) <- 6 + 5*X + alpha_Y[2]
p_Y[4] <- phi_Y[4]/sum(phi_Y)
log(phi_Y[4]) <- 7 + 6*X + alpha_Y[3]
alpha_Y ~ dmnorm.vcov(mu_alpha_Y, sigma_alpha_Y)
mu_alpha_Y[1] <- 0
sigma_alpha_Y[1, 1] <- 0.5
sigma_alpha_Y[1, 2] <- 0
sigma_alpha_Y[1, 3] <- 0
sigma_alpha_Y[1, 4] <- 0
mu_alpha_Y[2] <- 0
sigma_alpha_Y[2, 1] <- 0
sigma_alpha_Y[2, 2] <- 0.5
sigma_alpha_Y[2, 3] <- 0
sigma_alpha_Y[2, 4] <- 0
mu_alpha_Y[3] <- 0
sigma_alpha_Y[3, 1] <- 0
sigma_alpha_Y[3, 2] <- 0
sigma_alpha_Y[3, 3] <- 0.5
sigma_alpha_Y[3, 4] <- 0
mu_alpha_Y[4] <- 0
sigma_alpha_Y[4, 1] <- 0
sigma_alpha_Y[4, 2] <- 0
sigma_alpha_Y[4, 3] <- 0
sigma_alpha_Y[4, 4] <- 0.5")
  expect_output(categorical_bugsGroup(nodename = nodename,
                                      nodesCatIdx = nodesCatIdx,
                                      nodesintercepts = nodesintercepts,
                                      parentnames = parentnames,
                                      parentcoefs = parentcoefs,
                                      sigma = sigma,
                                      sigma_alpha = sigma_alpha),
                regexp = outstring,
                fixed = TRUE)
})

test_that("Poisson response Bugs model construction works", {
  if(.Platform$OS.type == "unix") {
    capture.output(
      # No parent nodes / no predictors
      expect_no_error(pois_bugs(nodename = "a",
                                parentnames = NULL,
                                nodesintercept = c(0.318077),
                                parentcoefs = NULL)),

      # One parent nodes / predictor
      expect_no_error(pois_bugs(nodename = "a",
                                parentnames = "b",
                                nodesintercept = c(0.318077),
                                parentcoefs = list("b"=c(b=0.3059395)))),

      # Multiple parent nodes / predictors
      expect_no_error(pois_bugs(nodename = "a",
                                parentnames = c("b", "c"),
                                nodesintercept = c(0.318077),
                                parentcoefs = list("b"=c(b=0.3059395),
                                                   "c"=c(c=0.5555)))),
      file = "/dev/null"
    )
  }

  outstring <- c(
    "a ~ dpois(lambda.a) # Poisson response
log(lambda.a) <- 0.318077 + 0.3059395*b + 0.5555*c # logistic regression")
  expect_output(pois_bugs(nodename = "a",
                          parentnames = c("b", "c"),
                          nodesintercept = c(0.318077),
                          parentcoefs = list("b"=c(b=0.3059395),
                                             "c"=c(c=0.5555))),
                regexp = outstring,
                fixed = TRUE)
})

test_that("Poisson response Bugs model construction with mixed-effects works", {
  muY <- 3
  sigma_alphaY <- 0.5
  betaX <- 4
  betaN <- 5
  if(.Platform$OS.type == "unix") {
    capture.output(
      # No parent nodes / no predictors
      expect_no_error(pois_bugsGroup(nodename = "Y",
                                     nodesintercept = muY,
                                     parentnames = NULL,
                                     parentcoefs = NULL,
                                     sigma_alpha = sigma_alphaY)),

      # One parent nodes / predictor
      expect_no_error(pois_bugsGroup(nodename = "Y",
                                     nodesintercept = muY,
                                     parentnames = c("X"),
                                     parentcoefs = list("X"=c(X=betaX)),
                                     sigma_alpha = sigma_alphaY)),

      # Multiple parent nodes / predictors
      expect_no_error(pois_bugsGroup(nodename = "Y",
                                     nodesintercept = muY,
                                     parentnames = c("X", "N"),
                                     parentcoefs = list("X"=c(X=betaX),
                                                        "N"=c(N=betaN)),
                                     sigma_alpha = sigma_alphaY)),

      # Multinomial parent nodes / predictors
      # expect_no_error(),
      file = "/dev/null"
    )
  }

  outstring <- c(
    "Y ~ dpois(lambda_Y)
log(lambda_Y) <- mu_Y + 4*X + 5*N + alpha_Y
mu_Y <- 3
alpha_Y ~ dnorm(mu_alpha_Y, tau_alpha_Y)
mu_alpha_Y <- 0
tau_alpha_Y <- inverse(sigma_alpha_Y)
sigma_alpha_Y <- 0.5")
  expect_output(pois_bugsGroup(nodename = "Y",
                               nodesintercept = muY,
                               parentnames = c("X", "N"),
                               parentcoefs = list("X"=c(X=betaX),
                                                  "N"=c(N=betaN)),
                               sigma_alpha = sigma_alphaY),
                regexp = outstring,
                fixed = TRUE)
})

test_that("makebugs() works", {
  mydists <- list(a="gaussian",
                  b="multinomial",
                  c="binomial",
                  d="poisson")
  mydag <- matrix(0, 4, 4, byrow = TRUE, dimnames = list(c("a", "b", "c", "d"), c("a", "b", "c", "d")))
  mydag[2,1] <- mydag[3,2] <- mydag[4,3] <- 1
  # plotAbn(mydag, data.dists = mydists)
  mycoefs <- list("a"=matrix(-6.883383e-17, byrow = T, dimnames = list(NULL, "a|intercept")),
                  "b"=matrix(c(2.18865, 3.133928, 3.138531, 1.686432, 3.134161, 5.052104),
                             nrow= 1, byrow = T, dimnames = list(c(NULL),
                                                                 c("b|intercept.2", "b|intercept.3", "b|intercept.4", "a.2", "a.3", "a.4"))),
                  "c"=matrix(c(1.11, 2.22, 3.33, 4.44, 5.55),
                             nrow= 1, byrow = T, dimnames = list(c(NULL),
                                                                 c("c|intercept", "b1", "b2", "b3", "b4"))),
                  "d"=matrix(c(3.33, 4.44),
                             nrow= 1, byrow = T, dimnames = list(c(NULL),
                                                                 c("d|intercept", "c"))))
  mymse <- c("a"=0,"b"=1,"c"=2,"d"=3)
  if(.Platform$OS.type == "unix") {
    capture.output(
      expect_no_error({
        makebugs(dag = mydag,
                 data.dists = mydists,
                 coefs = mycoefs,
                 std = mymse)}),
      file = "/dev/null")
  }


  outstring <- c(
    "model{
a ~ dnorm(mu.a, precision.a) # Gaussian response
mu.a <- -6.883383e-17 # Linear regression
precision.a <- inverse(0) # precision tau = 1/standard_dev
b ~ dcat(p.b) # Categorical response
p.b[1] <- phi.b[1]/sum(phi.b) # soft-max
log(phi.b[1]) <- 0 # Reference category
p.b[2] <- phi.b[2]/sum(phi.b) # soft-max
log(phi.b[2]) <- 2.18865 + 1.686432*a
p.b[3] <- phi.b[3]/sum(phi.b) # soft-max
log(phi.b[3]) <- 3.133928 + 3.134161*a
p.b[4] <- phi.b[4]/sum(phi.b) # soft-max
log(phi.b[4]) <- 3.138531 + 5.052104*a
c ~ dbern(p.c) # Bernoulli response
logit(p.c) <- 1.11 + 0.0242836145724376*b + 0.0736851897251164*b + 0.223587273987992*b + 0.678443921714454*b # logistic regression
d ~ dpois(lambda.d) # Poisson response
log(lambda.d) <- 3.33 + 4.44*c # logistic regression
}"
  )
  expect_output({
    makebugs(dag = mydag,
             data.dists = mydists,
             coefs = mycoefs,
             std = mymse)
  },
  regexp = outstring,
  fixed = TRUE,
  width = 150)
})

test_that("makebugsGroup() works",{
  # load("tests/testthat/testdata/makebugsGauss_data.Rdata")
  load("testdata/makebugsGauss_data.Rdata")

  expectedOut <- c(
    "model{
Outdoor ~ dbern(p_Outdoor)
logit(p_Outdoor) <- mu_Outdoor + alpha_Outdoor
mu_Outdoor <- -0.1231088
alpha_Outdoor ~ dnorm(mu_alpha_Outdoor, tau_alpha_Outdoor)
mu_alpha_Outdoor <- 0
tau_alpha_Outdoor <- inverse(sigma_alpha_Outdoor)
sigma_alpha_Outdoor <- 1.045681
Sex ~ dcat(p_Sex)
p_Sex[1] <- phi_Sex[1]/sum(phi_Sex)
log(phi_Sex[1]) <- 0 + alpha_Sex[1]
p_Sex[2] <- phi_Sex[2]/sum(phi_Sex)
log(phi_Sex[2]) <- 0.9706316 + 1.28692758032913*Age + alpha_Sex[1]
p_Sex[3] <- phi_Sex[3]/sum(phi_Sex)
log(phi_Sex[3]) <- -0.05845993 + 0.28240627414296*Age + alpha_Sex[2]
p_Sex[4] <- phi_Sex[4]/sum(phi_Sex)
log(phi_Sex[4]) <- 0.3357361 + 1.57473104113221*Age + alpha_Sex[3]
alpha_Sex ~ dmnorm.vcov(mu_alpha_Sex, sigma_alpha_Sex)
mu_alpha_Sex[1] <- 0
sigma_alpha_Sex[1, 1] <- 5.933619e-05
sigma_alpha_Sex[1, 2] <- 1e-50
sigma_alpha_Sex[1, 3] <- 4.32787e-08
mu_alpha_Sex[2] <- 0
sigma_alpha_Sex[2, 1] <- 1e-50
sigma_alpha_Sex[2, 2] <- 5.834134e-05
sigma_alpha_Sex[2, 3] <- 1e-50
mu_alpha_Sex[3] <- 0
sigma_alpha_Sex[3, 1] <- 4.32787e-08
sigma_alpha_Sex[3, 2] <- 1e-50
sigma_alpha_Sex[3, 3] <- 6.37651e-05
GroupSize ~ dpois(lambda_GroupSize)
log(lambda_GroupSize) <- mu_GroupSize + -0.282436260957026*Outdoor + alpha_GroupSize
mu_GroupSize <- 1.251637
alpha_GroupSize ~ dnorm(mu_alpha_GroupSize, tau_alpha_GroupSize)
mu_alpha_GroupSize <- 0
tau_alpha_GroupSize <- inverse(sigma_alpha_GroupSize)
sigma_alpha_GroupSize <- 1e-50
Age <- mu_Age + -0.00454841753065188*GroupSize + alpha_Age + e_Age
mu_Age <- -0.04832546
e_Age ~ dnorm(mu_e_Age, tau_Age)
mu_e_Age <- 0
tau_Age <- inverse(sigma_Age)
sigma_Age <- 0.9912502
alpha_Age ~ dnorm(mu_alpha_Age, tau_alpha_Age)
mu_alpha_Age <- 0
tau_alpha_Age <- inverse(sigma_alpha_Age)
sigma_alpha_Age <- 0.22154
}")

  expect_output({
    makebugsGroup(dag = dag,
                  data.dists = data.dists,
                  stderrors = mse,
                  group.var = group.var,
                  mu = mu,
                  betas = betas,
                  sigma = sigm,
                  sigma_alpha = sigm_alpha)},
    regexp = expectedOut,
    fixed = TRUE,
    width = 150)
})

test_that("simulateAbn() catches wrong arguments", {
  # Make a proper abnFit object
  ## without group.var
  if(.Platform$OS.type == "unix") {
    capture.output({
      df <- FCV[, c(12:15)]
      mydists <- list(Outdoor="binomial",
                      Sex="multinomial",
                      GroupSize="poisson",
                      Age="gaussian")

      ## buildScoreCache -> mostProbable() -> fitAbn()
      suppressWarnings({
        mycache.mle <- buildScoreCache(data.df = df, data.dists = mydists, method = "mle",
                                       adj.vars = NULL, cor.vars = NULL,
                                       dag.banned = NULL, dag.retained = NULL,
                                       max.parents = 1,
                                       which.nodes = NULL, defn.res = NULL)
      }) # ignore non-convergence warnings
      expect_no_error({
        mp.dag.mle <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
      })
      expect_no_error({
        myres.mle <- fitAbn(object = mp.dag.mle, method = "mle")
      })
    },
    file = "/dev/null"
    )

    # The actual tests
    capture.output(
      expect_no_error({
        simulateAbn(object = myres.mle,
                    run.simulation = TRUE,
                    bugsfile = NULL,
                    verbose = FALSE)
      }),
      file = "/dev/null")
    expect_error({
      simulateAbn(object = myres.mle,
                  run.simulation = "TRUE",
                  bugsfile = NULL,
                  verbose = FALSE)
    })
    expect_error({
      simulateAbn(object = myres.mle,
                  run.simulation = TRUE,
                  bugsfile = 1,
                  verbose = FALSE)
    })
    expect_error({
      simulateAbn(object = myres.mle,
                  run.simulation = TRUE,
                  bugsfile = NULL,
                  verbose = "FALSE")
    })
    expect_error({
      simulateAbn(object = unclass(myres.mle),
                  run.simulation = TRUE,
                  bugsfile = NULL,
                  verbose = FALSE)
    })
    expect_error({
      simulateAbn(object = NULL,
                  run.simulation = TRUE,
                  bugsfile = NULL,
                  verbose = FALSE)
    })
  }
})

test_that("simulateAbn() simulation works with method 'mle'", {
  # Make a proper abnFit object
  ## without group.var
  if(.Platform$OS.type == "unix") {
    capture.output({
      df <- FCV[, c(12:15)]
      mydists <- list(Outdoor="binomial",
                      Sex="multinomial",
                      GroupSize="poisson",
                      Age="gaussian")

      ## buildScoreCache -> mostProbable() -> fitAbn()
      suppressWarnings({
        mycache.mle <- buildScoreCache(data.df = df, data.dists = mydists, method = "mle",
                                       adj.vars = NULL, cor.vars = NULL,
                                       dag.banned = NULL, dag.retained = NULL,
                                       max.parents = 1,
                                       which.nodes = NULL, defn.res = NULL)
      }) # ignore non-convergence warnings
      expect_no_error({
        mp.dag.mle <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
      })
      expect_no_error({
        myres.mle <- fitAbn(object = mp.dag.mle, method = "mle", centre = FALSE)
      })
    },
    file = "/dev/null"
    )

    expect_no_error({
      myres.sim <- simulateAbn(object = myres.mle,
                               run.simulation = TRUE,
                               bugsfile = NULL,
                               n.chains = 10L,
                               n.adapt = 1000L,
                               n.thin = 100L,
                               n.iter = 10000L,
                               seed = 42L,
                               verbose = FALSE)
    })

    act <- as.numeric(round(prop.table(table(myres.sim$Outdoor)), 2))
    expected <- as.numeric(round(prop.table(table(df$Outdoor)), 2))
    expect_equal(act[order(act)],
                 expected[order(expected)],
                 tolerance = 0.05)
    act <- as.numeric(round(prop.table(table(myres.sim$Sex)), 2))
    expected <- as.numeric(round(prop.table(table(df$Sex)), 2))
    expect_equal(act[order(act)],
                 expected[order(expected)],
                 tolerance = 0.05)

    ## with group.var
    suppressWarnings({
      suppressMessages({
        capture.output({
          df <- FCV[, c(11:15)]
          mydists <- list(Pedigree="binomial",
                          Outdoor="binomial",
                          Sex="multinomial",
                          GroupSize="poisson",
                          Age="gaussian")
          mydists <- mydists[-1] # remove grouping variable from distribution list
          retaindag <- matrix(0, nrow = length(mydists), ncol = length(mydists), dimnames = list(names(mydists), names(mydists)))
          retaindag[3, 1] <- retaindag[4,3] <- retaindag[2,4] <- 1
          mycache.mle.grp <- buildScoreCache(data.df = df, data.dists = mydists, method = "mle",
                                             group.var = "Pedigree", adj.vars = NULL, cor.vars = NULL,
                                             dag.banned = NULL, dag.retained = retaindag,
                                             max.parents = 3,
                                             which.nodes = NULL, defn.res = NULL)
          mp.dag.mle.grp <- mostProbable(score.cache = mycache.mle.grp, verbose = FALSE)
          expect_no_error({
            myres.mle.grp <- fitAbn(object = mp.dag.mle.grp, method = "mle", group.var = "Pedigree")
          })
          expect_no_error({
            mysim.grp <- simulateAbn(object = myres.mle.grp,
                                     run.simulation = TRUE,
                                     bugsfile = NULL,
                                     verbose = FALSE,
                                     debug = FALSE)
          })
        },
        file = "/dev/null"
        )
      })
    })
    act <- as.numeric(round(prop.table(table(mysim.grp$Outdoor)), 2))
    expected <- as.numeric(round(prop.table(table(df$Outdoor)), 2))
    expect_equal(act[order(act)],
                 expected[order(expected)],
                 tolerance = 0.1) # quite high tolerance because the actual data has no grouping...

    # Correct number of categories simulated for categorical/multinomial variables?
    act <- as.numeric(round(prop.table(table(mysim.grp$Sex)), 2))
    expected <- as.numeric(round(prop.table(table(df$Sex)), 2))
    expect_equal(length(act), length(expected))
  }
})

test_that("simulateAbn() simulation works with method 'bayes'", {
  # Make a proper abnFit object
  ## without group.var
  if(.Platform$OS.type == "unix") {
    capture.output({
      df <- FCV[, c(12, 14:15)]
      mydists <- list(Outdoor="binomial",
                      # Sex="multinomial",
                      GroupSize="poisson",
                      Age="gaussian")

      ## buildScoreCache -> mostProbable() -> fitAbn()
      suppressWarnings({
        mycache.bayes <- buildScoreCache(data.df = df, data.dists = mydists, method = "bayes",
                                       adj.vars = NULL, cor.vars = NULL,
                                       dag.banned = NULL, dag.retained = NULL,
                                       max.parents = 1,
                                       which.nodes = NULL, defn.res = NULL)
      }) # ignore non-convergence warnings
      expect_no_error({
        mp.dag.bayes <- mostProbable(score.cache = mycache.bayes, verbose = FALSE)
      })
      expect_no_error({
        myres.bayes <- fitAbn(object = mp.dag.bayes, method = "bayes", centre = FALSE)
      })
    },
    file = "/dev/null"
    )

    expect_no_error({
      myres.sim <- simulateAbn(object = myres.bayes,
                               run.simulation = TRUE,
                               bugsfile = NULL,
                               n.chains = 10L,
                               n.adapt = 1000L,
                               n.thin = 100L,
                               n.iter = 10000L,
                               seed = 42L,
                               verbose = FALSE)
    })

    act <- as.numeric(round(prop.table(table(myres.sim$Outdoor)), 2))
    expected <- as.numeric(round(prop.table(table(df$Outdoor)), 2))
    expect_equal(act[order(act)],
                 expected[order(expected)],
                 tolerance = 0.05)
    act <- as.numeric(round(prop.table(table(myres.sim$Sex)), 2))
    expected <- as.numeric(round(prop.table(table(df$Sex)), 2))
    expect_equal(act[order(act)],
                 expected[order(expected)],
                 tolerance = 0.05)

    skip("simulateAbn() with method 'bayes' is not tested with 'group.var'.")
    ## with group.var
    suppressWarnings({
      suppressMessages({
        capture.output({
          df <- FCV[, c(11:15)]
          mydists <- list(Pedigree="binomial",
                          Outdoor="binomial",
                          Sex="multinomial",
                          GroupSize="poisson",
                          Age="gaussian")
          mydists <- mydists[-1] # remove grouping variable from distribution list
          retaindag <- matrix(0, nrow = length(mydists), ncol = length(mydists), dimnames = list(names(mydists), names(mydists)))
          retaindag[3, 1] <- retaindag[4,3] <- retaindag[2,4] <- 1
          mycache.mle.grp <- buildScoreCache(data.df = df, data.dists = mydists, method = "bayes",
                                             group.var = "Pedigree", adj.vars = NULL, cor.vars = NULL,
                                             dag.banned = NULL, dag.retained = retaindag,
                                             max.parents = 3,
                                             which.nodes = NULL, defn.res = NULL)
          mp.dag.mle.grp <- mostProbable(score.cache = mycache.mle.grp, verbose = FALSE)
          expect_no_error({
            myres.mle.grp <- fitAbn(object = mp.dag.mle.grp, method = "bayes", group.var = "Pedigree")
          })
          expect_no_error({
            mysim.grp <- simulateAbn(object = myres.mle.grp,
                                     run.simulation = TRUE,
                                     bugsfile = NULL,
                                     verbose = FALSE,
                                     debug = FALSE)
          })
        },
        file = "/dev/null"
        )
      })
    })
    act <- as.numeric(round(prop.table(table(mysim.grp$Outdoor)), 2))
    expected <- as.numeric(round(prop.table(table(df$Outdoor)), 2))
    expect_equal(act[order(act)],
                 expected[order(expected)],
                 tolerance = 0.1) # quite high tolerance because the actual data has no grouping...

    # Correct number of categories simulated for categorical/multinomial variables?
    act <- as.numeric(round(prop.table(table(mysim.grp$Sex)), 2))
    expected <- as.numeric(round(prop.table(table(df$Sex)), 2))
    expect_equal(length(act), length(expected))
  }
})

test_that("simulateAbn() works with grouping in real data.",{
  suppressMessages({
    suppressWarnings({
      if(.Platform$OS.type == "unix") {
        capture.output({
          df <- adg
          df[,1:5] <- lapply(df[,1:5], factor)
          df[,9] <- factor(df[,9])
          mydists <- list(AR = "binomial",
                          pneumS = "binomial",
                          female = "binomial",
                          livdam = "binomial",
                          eggs = "binomial",
                          wormCount = "poisson",
                          age = "gaussian",
                          adg = "gaussian")
          # farm = "multinomial")
          #ban/retain matrices
          myretain <- matrix(0, length(mydists), length(mydists))
          colnames(myretain) <- rownames(myretain) <- names(mydists)

          mybanned <- matrix(0, length(mydists), length(mydists))
          colnames(mybanned) <- rownames(mybanned) <- names(mydists)

          mybanned[3,-3] <- 1

          ###
          # max parent search
          ###
          max.par <- length(mydists)-1 # maximal possible parent values
          all.fits <- list()
          for (i in 1:max.par) {
            mycache <- buildScoreCache(data.df = df,
                                       data.dists = mydists,
                                       dag.banned = mybanned,
                                       dag.retained = myretain,
                                       max.parents = i,
                                       group.var = "farm",
                                       method = "mle")
            mydag <- mostProbable(score.cache = mycache)
            fabn <- fitAbn(object = mydag, method = "mle")
            cat(paste("network score for", i, "parents =", fabn$mlik, "\n\n"))
            all.fits[i] <- list(fabn)
          }
          allmlik <- lapply(X = all.fits, FUN = function(x){sum(x$mlik)})
          # plot(unlist(allmlik), xlab = "number of max. parents", ylab = "mlik")
          first.maxpar <- min(which(unlist(allmlik) == max(unlist(allmlik))))

          ###
          # extract best fit
          ###
          best.fit <- all.fits[[first.maxpar]]
          # plotAbn(best.fit)

          ###
          # simulate new data based on estimated parameters from GLMM
          ###
          mysim <- simulateAbn(best.fit,
                               verbose = FALSE)

          act <- as.numeric(round(prop.table(table(mysim$livdam)), 2))
          expected <- as.numeric(round(prop.table(table(df$livdam)), 2))
          expect_equal(act[order(act)],
                       expected[order(expected)],
                       tolerance = 0.001)

          act <- as.numeric(round(prop.table(table(mysim$eggs)), 2))
          expected <- as.numeric(round(prop.table(table(df$eggs)), 2))
          expect_equal(act[order(act)],
                       expected[order(expected)],
                       tolerance = 0.05)

          act <- as.numeric(round(prop.table(table(mysim$pneumS)), 2))
          expected <- as.numeric(round(prop.table(table(df$pneumS)), 2))
          expect_equal(act[order(act)],
                       expected[order(expected)],
                       tolerance = 0.07) # quite high tolerance because the actual data has no grouping...
        },
        file = "/dev/null"
        )
      }
    })
  })
})
