
# abn - Additive Bayesian Networks

<!-- badges: start -->
[![pipeline status](https://git.math.uzh.ch/mdeluc/devel-abn/badges/master/pipeline.svg?key_text=R-CMD-check&key_width=90)](https://git.math.uzh.ch/mdeluc/devel-abn/-/commits/master)
[![coverage report](https://git.math.uzh.ch/mdeluc/devel-abn/badges/master/coverage.svg)](https://git.math.uzh.ch/mdeluc/devel-abn/-/commits/master)
[![Latest Release](https://git.math.uzh.ch/mdeluc/devel-abn/-/badges/release.svg)](https://git.math.uzh.ch/mdeluc/devel-abn/-/releases)
![cran](https://www.r-pkg.org/badges/version-ago/abn) 
![downloads](https://cranlogs.r-pkg.org/badges/grand-total/abn) 
<!-- badges: end -->

Bayesian network analysis is a form of probabilistic graphical models which derives from empirical data a directed acyclic graph, DAG, describing the dependency structure between random variables.

An additive Bayesian network model consists of a form of a DAG where each node comprises a generalized linear model, GLM. Additive Bayesian network models are equivalent to Bayesian multivariate regression using graphical modelling, they generalises the usual multivariable regression, GLM, to multiple dependent variables.

'abn' provides routines to help determine optimal Bayesian network models for a given data set, where these models are used to identify statistical dependencies in messy, complex data. The additive formulation of these models is equivalent to multivariate generalised linear modelling (including mixed models with iid random effects).

The usual term to describe this model selection process is structure discovery.

The core functionality is concerned with model selection - determining the most robust empirical model of data from interdependent variables. Laplace approximations are used to estimate goodness of fit metrics and model parameters, and wrappers are also included to the INLA package which can be obtained from <https://www.r-inla.org>.

The computing library JAGS <https://mcmc-jags.sourceforge.io> is used to simulate 'abn'-like data.

A comprehensive set of documented case studies, numerical accuracy/quality assurance exercises, and additional documentation are available from the 'abn' website <http://r-bayesian-networks.org>.

## Installation

You can install the latest release from CRAN with 
``` r
install.packages("abn")
```

The recent development version can be installed from:
``` r
remotes::install_gitlab("https://git.math.uzh.ch/mdeluc/abn.git")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(abn)
## Built-in dataset with a subset of cols
mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4", "g2")]

## setup distribution list for each node
mydists <- list(b1="binomial", 
                b2="binomial", 
                b3="binomial",
                g1="gaussian",
                b4="binomial",
                p2="poisson",
                p4="poisson",
                g2="gaussian")

# Structural constraints
# ban arc from b2 to b1
# always retain arc from g2 to g1

## now build the cache of pre-computed scores accordingly to the structural constraints
mycache.mle <- buildScoreCache(data.df=mydat, 
                         data.dists=mydists,
                         dag.banned= ~b1|b2, 
                         dag.retained= ~g1|g2, 
                         max.parents=3,
                         method = "mle")

## Find most probable DAG
mp.dag.mle <- mostProbable(mycache.mle)

## Now fit the model to calculate its goodness-of-fit
myres <- fitAbn(object=mp.dag.mle,
                method = "mle")

## Plot the DAG
plotAbn(myres)

## Log-marginal likelihood goodness-of-fit for complete DAG
print(myres$mlik)
```
