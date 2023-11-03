
# ABN: Additive Bayesian Networks

<!-- badges: start -->
[![pipeline status](https://git.math.uzh.ch/mdeluc/abn/badges/main/pipeline.svg?key_text=R-CMD-check&key_width=90)](https://git.math.uzh.ch/mdeluc/abn/-/commits/main)
[![coverage report](https://git.math.uzh.ch/mdeluc/abn/badges/main/coverage.svg)](https://git.math.uzh.ch/mdeluc/abn/-/commits/main)
[![Latest Release](https://git.math.uzh.ch/mdeluc/abn/-/badges/release.svg)](https://git.math.uzh.ch/mdeluc/abn/-/commits/main)
![cran](https://www.r-pkg.org/badges/version-ago/abn) 
![downloads](https://cranlogs.r-pkg.org/badges/grand-total/abn) 
![LICENCE](https://img.shields.io/cran/l/abn)
<!-- badges: end -->

The R package `abn` is a tool for Bayesian network analysis, a form of probabilistic graphical models. It derives a directed acyclic graph (DAG) from empirical data, that describes the dependency structure between random variables. The package provides routines for structure learning and parameter estimation of additive Bayesian network models.  

## Installation

You can install the latest release from CRAN as follows:

``` r
install.packages("abn")
```

The recent development version can be installed from:
``` r
remotes::install_gitlab("https://git.math.uzh.ch/mdeluc/abn.git")
```

## What is an additive Bayesian network?

Additive Bayesian network (ABN) models are a type of statistical model that uses the principles of Bayesian statistics and graph theory. They provide a framework for representing data with multiple variables, known as multivariate data.

ABN models can be seen as a graphical representation of (Bayesian) multivariate regression. 
This form of statistical analysis enables the prediction of multiple outcomes from a given set of predictors, while simultaneously accounting for the relationships between these outcomes.

In other words, additive Bayesian network models extend the concept of Generalized Linear Models (GLMs), which are typically used for predicting a single outcome, to scenarios where there are multiple dependent variables. This makes them a powerful tool for understanding complex, multivariate datasets.

## Features

The R package `abn` provides routines to help determine optimal additive Bayesian network models for a given data set. 
The core functionality is concerned with model selection - determining the most likely model of data from interdependent variables. 
Expert knowledge can be incorporated into the model selection process by specifying structural constraints, such as which arcs are banned or retained.

The general workflow with `abn` follows a three step process:
1. **Determine the model search space**: The function `buildScoreCache()` builds a cache of pre-computed scores for each possible DAG.
This requires the specification of the data set, the data types of the variables, and the structural constraints of the model (e.g. which arcs are banned or retained and the maximum number of parents per node).

2. **Structure learning**: `abn` offers different structure learning algorithms:
    - The exact structure learning algorithm from [Koivisto and Sood (2004)](https://www.jmlr.org/papers/volume5/koivisto04a/koivisto04a.pdf) is implemented in `C` and can be called with the function `mostProbable()`, which finds the most probable DAG for a given data set.
    - A set of heuristic search algorithms is available through the function `searchHeuristic()`. It includes the hill-climber, tabu search, and simulated annealing algorithms that are implemented in `R`.
    - `searchHillClimber()` searches for high scoring DAGs using a random re-start greedy hill-climber heuristic search and is implemented in `C`. It slightly deviates from the method initially presented by [Heckerman et al. 1995](https://doi.org/10.1023/A:1022623210503) (for details consult the respective help page `?abn::searchHillClimber()`).

3. **Parameter estimation**: With the function `fitAbn()` the parameters of the model are estimated based on the DAG from the previous step.

`abn` allows for two different model formulations, specified with the argument `method`:
- `method = "mle"` fits a model under the frequentist paradigm using information-theoretic criteria to select the best model.
- `method = "bayes"` estimates the posterior distribution of the model parameters based on two Laplace approximation methods, which is a method for Bayesian inference and an alternative to Markov Chain Monte Carlo (MCMC). A standard Laplace approximation is implemented in the `abn` source code but switches in specific cases (see help page `?abn::fitAbn()`) to the Integrated Nested Laplace Approximation from the [INLA package](https://www.r-inla.org) requiring the installation thereof.

To generate new observations from a fitted ABN model, the function `simulateAbn()` that simulate data based on the DAG and the estimated parameters from the previous step. `simulateAbn()` is available for both `method = "mle"` and `method = "bayes"` and requires the installation of the [JAGS package](https://mcmc-jags.sourceforge.io). 

### Supported Data types

The `abn` package supports the following distributions for the variables in the network:
- Gaussian distribution for continuous variables.
- Binomial distribution for binary variables.
- Poisson distribution for variables with count data.
- Multinomial distribution for categorical variables (only available with `method = "mle"`).

Unlike other packages, `abn` does not restrict the combination of parent-child distributions.

### Multilevel Models for Grouped Data Structures

The analysis of "hierarchical" or "grouped" data, in which observations are nested within higher-level units, requires statistical models with parameters that vary across groups (e.g. mixed-effect models).

`abn` allows to control for one layer clustering, where observations are grouped into a single layer of clusters which are themself assumed to be independent, but observations within the clusters may be correlated (e.g. students nested within schools, measurements over time for each patient, etc).
The argument `group.var` specifies the discrete variable that defines the group structure. The model is then fitted separately for each group and the results are combined. 

For example, studying student test scores across different schools, a varying intercept model would allow for the possibility that average test scores (the intercept) might be higher in one school than another due to factors specific to each school. This can be modeled in `abn` with setting the argument `group.var` to the variable containing the school names. The model is then fitted as varying intercept model, where the intercept is allowed to vary across schools but the slope is assumed to be the same for all schools.

Under the frequentist paradigm (`method = "mle"`), `abn` relies on the `lme4` package to fit generalized linear mixed models (GLMMs) for Binomial, Poisson, and Gaussian distributed variables. For multinomial distributed variables, `abn` fits a multinomial baseline category logit model with random effects using the `mclogit` package. Currently only one layer clustering is supported (e.g. for `method = "mle"` this corresponds to a random intercept model)`.

With a Bayesian approach (`method = "bayes"`), `abn` relies on its own implementation of the Laplace approximation and the package `INLA` to fit a single level hierarchical model for Binomial, Poisson, and Gaussian distributed variables. Multinomial distributed variables in general (see Section [Data Types](#data-types)) are not yet implemented with `method = "bayes"`.

### Applications

A comprehensive set of documented case studies, numerical accuracy/quality assurance exercises, and additional documentation are available from the [`abn` website](http://r-bayesian-networks.org).

## Examples

- [Example 1: Basic usage](#example-1-basic-usage)
- [Example 2: Restrict model search space](#example-2-restrict-model-search-space)
- [Example 3: Grouped Data Structures](#example-3-grouped-data-structures)
- [Example 4: Using INLA vs internal Laplace approximation](#example-4-using-inla-vs-internal-laplace-approximation)

### Example 1: Basic Usage

This is a basic example which shows the basic workflow:

``` r
library(abn)

# Built-in toy dataset with two Gaussian variables G1 and G2, two Binomial variables B1 and B2, and one multinomial variable C
str(g2b2c_data)

# Define the distributions of the variables
dists <- list(G1 = "gaussian",
              B1 = "binomial",
              B2 = "binomial",
              C = "multinomial",
              G2 = "gaussian")


# Build the score cache
cacheMLE <- buildScoreCache(data.df = g2b2c_data,
                         data.dists = dists,
                         method = "mle",
                         max.parents = 2)

# Find the most probable DAG
dagMP <- mostProbable(score.cache = cacheMLE)

# Print the most probable DAG
dagMP$dag

# Plot the most probable DAG
plot(dagMP)

# Fit the most probable DAG
myfit <- fitAbn(object = dagMP,
              method = "mle")

# Print the fitted DAG
print(myfit)

```

### Example 2: Restrict Model Search Space

Based on [example 1](#example-1-basic-usage) we may know that the arc G1->G2 is not possible and that the arc from C -> G2 must be present.
This "expert knowledge" can be included in the model by banning the arc from b2 to b1 and retaining the arc from g2 to g1.
The ratain and ban matrices are specified as an adjacency matrix of 0 and 1 entries, where 1 indicates that the arc is banned or retained, respectively. 
Row and column names must match the variable names in the data set. 
The corresponding column is a parent of the variable in the row.
Each column represent the parents and the rows the children.

Further we may want to restrict the maximum number of parents per node to 2.

``` r
# Ban the edge G1 -> G2
banmat <- matrix(0, nrow = 5, ncol = 5, dimnames = list(names(dists), names(dists)))
banmat[1, 5] <- 1

# retain always the edge C -> G2
retainmat <- matrix(0, nrow = 5, ncol = 5, dimnames = list(names(dists), names(dists)))
retainmat[5, 4] <- 1

# Limit the maximum number of parents to 2
max.par <- 2

# Build the score cache
cacheMLE_small <- buildScoreCache(data.df = g2b2c_data,
                            data.dists = dists,
                            method = "mle",
                            dag.banned = banmat,
                            dag.retained = retainmat,
                            max.parents = max.par)
print("Without restrictions from example 1: ", nrow(cacheMLE$node.defn))
print("With restrictions as in example 2: ", nrow(cacheMLE_small$node.defn))
```

### Example 3: Grouped Data Structures

Depending on the data structure, we may want to control for one layer clustering, where observations are grouped into a single layer of clusters which are themself assumed to be independent, but observations within the clusters may be correlated (e.g. students nested within schools, measurements over time for each patient, etc).

Currently, `abn` supports only one layer clustering. 

``` r

# Built-in toy data set
str(g2pbcgrp)

# Define the distributions of the variables
dists <- list(G1 = "gaussian",
              P = "poisson",
              B = "binomial",
              C = "multinomial",
              G2 = "gaussian") # group is not among the list of variable distributions

# Ban arcs such that C has only B and P as parents
ban.mat <- matrix(0, nrow = 5, ncol = 5, dimnames = list(names(dists), names(dists)))
ban.mat[4, 1] <- 1
ban.mat[4, 4] <- 1
ban.mat[4, 5] <- 1

# Build the score cache
cache <- buildScoreCache(data.df = g2pbcgrp,
                         data.dists = dists,
                         group.var = "group",
                         dag.banned = ban.mat,
                         method = "mle",
                         max.parents = 2)

# Find the most probable DAG
dag <- mostProbable(score.cache = cache)

# Plot the most probable DAG
plot(dag)

# Fit the most probable DAG
fit <- fitAbn(object = dag,
              method = "mle")

# Plot the fitted DAG
plot(fit)

# Print the fitted DAG
print(fit)
```


### Example 4: Using INLA vs internal Laplace approximation

Under a Bayesian approach, `abn` automatically switches to the Integrated Nested Laplace Approximation from the [INLA package](https://www.r-inla.org) if the internal Laplace approximation fails to converge. 
However, we can also force the use of INLA by setting the argument `control=list(max.mode.error=100)`.
The following example shows that the results are very similar and it also shows how to constrain arcs as formula objects and how to specify different parent limits for each node separately.

``` r
library(abn)

# Subset of the build-in dataset, see  ?ex0.dag.data
mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")] ## take a subset of cols

# setup distribution list for each node
mydists <- list(b1="binomial", b2="binomial", g1="gaussian",
                g2="gaussian", b3="binomial", g3="gaussian")

# Structural constraints
## ban arc from b2 to b1
## always retain arc from g2 to g1
## parent limits - can be specified for each node separately
max.par <- list("b1"=2, "b2"=2, "g1"=2, "g2"=2, "b3"=2, "g3"=2)

# now build the cache of pre-computed scores accordingly to the structural constraints
res.c <- buildScoreCache(data.df=mydat, data.dists=mydists,
                         dag.banned= ~b1|b2, 
                         dag.retained= ~g1|g2, 
                         max.parents=max.par)


# repeat but using R-INLA. The mlik's should be virtually identical.
if(requireNamespace("INLA", quietly = TRUE)){
  res.inla <- buildScoreCache(data.df=mydat, data.dists=mydists,
                              dag.banned= ~b1|b2, # ban arc from b2 to b1
                              dag.retained= ~g1|g2, # always retain arc from g2 to g1
                              max.parents=max.par,
                              control=list(max.mode.error=100)) # force using of INLA
  
  ## comparison - very similar
  difference <- res.c$mlik - res.inla$mlik
  summary(difference)
}
```

## Contributing

We welcome contributions to the `abn` package directly via pull requests.
Please motivate your contribution first in an issue, so we can discuss the best way to integrate it into the package.

## Citation

If you use `abn` in your research, please cite it as follows:

``` r
> citation("abn")
To cite abn in publications use:

  Kratzer G, Lewis F, Comin A, Pittavino M, Furrer R (2023). “Additive Bayesian Network Modeling with the R Package abn.” _Journal of Statistical Software_, *105*(8), 1-41.
  doi:10.18637/jss.v105.i08 <https://doi.org/10.18637/jss.v105.i08>.

To cite an example of a typical ABN analysis in publications use:

  Kratzer, G., Lewis, F.I., Willi, B., Meli, M.L., Boretti, F.S., Hofmann-Lehmann, R., Torgerson, P., Furrer, R. and Hartnack, S. (2020). Bayesian Network Modeling Applied to
  Feline Calicivirus Infection Among Cats in Switzerland. Frontiers in Veterinary Science, 7, 73

To cite package of the R package 'abn' in publications use:

  Furrer, R., Kratzer, G. and Lewis, F.I. (2023). abn: Modelling Multivariate Data with Additive Bayesian Networks. R package version 2.7-2.
  https://CRAN.R-project.org/package=abn
```

## License

The `abn` package is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Code of Conduct

Please note that the `abn` project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
