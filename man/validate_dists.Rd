% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{validate_dists}
\alias{validate_dists}
\title{Check for valid distribution}
\usage{
validate_dists(data.dists, returnDists = TRUE, ...)
}
\arguments{
\item{data.dists}{list of variable distributions.}

\item{returnDists}{if TRUE (default) returns the same list as provided.}

\item{...}{additional arguments.}
}
\value{
either TRUE/FALSE or list of variable distributions as provided as input.
}
\description{
The distribution names must match \verb{inla() family=''}.
Similar to \code{get.var.types()}, mainly different in output.
}
\keyword{internal}
