% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulateabn.R
\name{pois_bugs}
\alias{pois_bugs}
\alias{pois_bugsGroup}
\title{Bugs code for Poisson response}
\usage{
pois_bugs(nodename, nodesintercept, parentnames, parentcoefs)

pois_bugsGroup(nodename, nodesintercept, parentnames, parentcoefs, sigma_alpha)
}
\arguments{
\item{nodename}{character string of response variable name.}

\item{nodesintercept}{overall mean of response. Parameter from fixed-effects intercept.}

\item{parentnames}{single character string (for one parent) or vector of characters (for multiple parent nodes) with parent node (predictor variables) names.}

\item{parentcoefs}{overall slope for each predictor (parent node) variable (fixed-effects).}

\item{sigma_alpha}{between-group variance. Parameter from random-effects intercept.}
}
\value{
Bugs model returned as stdout.
}
\description{
Bugs model for count response variable
\eqn{X \sim \mathcal{Pois}(\lambda)}.
}
\section{Functions}{
\itemize{
\item \code{pois_bugsGroup()}: Bugs code for Poisson response with varying intercept

}}
\examples{
pois_bugs(nodename = "a",
          parentnames = c("b", "c"),
          nodesintercept = c(0.318077),
          parentcoefs = list("b"=c(b=0.3059395),
                             "c"=c(c=0.5555)))
}
\seealso{
\link{makebugs} \link{simulateAbn}
}
\keyword{internal}
\keyword{utilities}
