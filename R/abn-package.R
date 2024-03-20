#' @title \code{abn} Package
#' @description
#' \code{abn} is a collection of functions for fitting, selecting/learning, analyzing, reporting Additive Bayesian Networks.
#' @section General overview:
#' What is \pkg{abn}:
#' Bayesian network modeling is a data analysis technique that is ideally suited to messy, highly correlated, and complex data sets.
#' This methodology is somewhat distinct from other forms of statistical modeling in that its focus is on structure discovery - determining an optimal graphical model that describes the inter-relationships in the underlying processes which generated the data.
#' It is a multivariate technique and can be used for one or many dependent variables.
#' This is a data-driven approach, as opposed to, rely only on subjective expert opinion to determine how variables of interest are inter-related (for example, structural equation modeling).
#'
#' The R package \pkg{abn} is designed to fit additive Bayesian models to observational data sets.
#' It contains routines to score Bayesian Networks based on Bayesian or information-theoretic formulation of generalized linear models.
#' It is equipped with exact search and greedy search algorithms to select the best network.
#' The Bayesian implementation supports random effects to control for one layer clustering.
#' It supports a possible mixture of continuous, discrete, and count data and inputs of prior knowledge at a structural level.
#'
#' The R package \pkg{abn} requires the R package \pkg{Rgraphviz} to work well.
#' It is store outside of CRAN; see \sQuote{Examples} for the code to install the last version.
#'
#' The web page \href{https://r-bayesian-networks.org/}{r-bayesian-networks.org} provides further case studies.
#' See also the files provided in the package directories \code{inst/bootstrapping_example} and \code{inst/old_vignette} for more details.
#'
#' @references
#' Kratzer, Gilles, Fraser Lewis, Arianna Comin, Marta Pittavino, and Reinhard Furrer. “Additive Bayesian Network Modeling with the R Package Abn.” Journal of Statistical Software 105 (January 28, 2023): 1–41. https://doi.org/10.18637/jss.v105.i08.
#'
#' Kratzer, G., Lewis, F.I., Comin, A., Pittavino, M. and Furrer, R. (2019). "Additive Bayesian Network Modelling with the R Package abn". arXiv preprint arXiv:1911.09006.
#'
#' Lewis, F. I., and Ward, M. P. (2013). "Improving epidemiologic data analyses through multivariate regression modeling". Emerging themes in epidemiology, 10(1), 4.
#'
#' Kratzer, G., Pittavino, M, Lewis, F. I., and Furrer, R., (2017). "abn: an R package for modelling multivariate data using additive Bayesian networks". R package version 2.2.  https://CRAN.R-project.org/package=abn
#'
#' @examples
#' ## Citations:
#' print(citation('abn'), bibtex=TRUE)
#'
#' ## Installing the R package Rgraphviz:
#' # if (!requireNamespace("BiocManager", quietly = TRUE))
#' #     install.packages("BiocManager")
#' # BiocManager::install("Rgraphviz")
#'
#' ## README.md in the directory `bootstrapping_example/`:
#' # edit(file=paste0( path.package('abn'),'/bootstrapping_example/README.md'))
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @returns nothing.
#' @useDynLib abn, .registration=TRUE
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' @title Documentation of C Functions
#' @description This is mainly to circumvent issues in R CMD check.
#' @returns nothing.
#' @name Cfunctions
#' @aliases buildcachematrix checkforcycles fit_single_node fitabn_marginals mostprobable_C searchhill
#' @keywords internal
NULL
