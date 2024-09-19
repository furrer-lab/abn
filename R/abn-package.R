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
#' Kratzer G, Lewis F, Comin A, Pittavino M, Furrer R (2023). “Additive Bayesian Network Modeling with the R Package abn.” Journal of Statistical Software, 105(8), 1-41. https://doi.org/10.18637/jss.v105.i08.
#'
#' Delucchi M, Liechti J, Spinner G, Furrer R (2024). “abn: Additive Bayesian Networks.” Journal of Open Source Software, Accepted for publication. R package version 3.1.3, https://joss.theoj.org/papers/1bbc43a2be86f5d3f831cedb5cf81812.
#'
#' Lewis, F. I., and Ward, M. P. (2013). "Improving epidemiologic data analyses through multivariate regression modeling". Emerging themes in epidemiology, 10(1), 4.
#'
#' Kratzer, G., Lewis, F.I., Willi, B., Meli, M.L., Boretti, F.S., Hofmann-Lehmann, R., Torgerson, P., Furrer, R. and Hartnack, S. (2020). "Bayesian Network Modeling Applied to Feline Calicivirus Infection Among Cats in Switzerland". Frontiers in Veterinary Science, 7, 73.
#'
#' Delucchi M, Furrer R, Kratzer G, Lewis F, Liechti J, Pittavino M, Cherneva K (2024). "abn: Modelling Multivariate Data with Additive Bayesian Networks". R package version 3.1.3, https://CRAN.R-project.org/package=abn.
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
