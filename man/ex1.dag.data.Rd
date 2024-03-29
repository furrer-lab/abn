% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ex1.dag.data}
\alias{ex1.dag.data}
\title{Synthetic validation data set for use with abn library examples}
\format{
A data frame, binary variables are factors.
The relevant formulas are given below (note these do not give parameter
estimates just the form of the relationships, like in glm(),
e.g. logit()=1+p1 means a logit link function and comprises of an
intercept term and a term involving p1).
\describe{
\item{b1}{binary, logit()=1 }
\item{p1}{poisson, log()=1 }
\item{g1}{gaussian, identity()=1 }
\item{b2}{binary, logit()=1}
\item{p2}{poisson, log()=1+b1+p1 }
\item{b3}{binary, logit()=1+b1+g1+b2 }
\item{g2}{gaussian, identify()=1+p1+g1+b2 }
\item{b4}{binary, logit()=1+g1+p2}
\item{b5}{binary, logit()=1+g1+g2 }
\item{g3}{gaussian, identity()=1+g1+b2 }
}
}
\usage{
ex1.dag.data
}
\description{
10000 observations simulated from a DAG with 10 variables from Poisson, Bernoulli and Gaussian distributions.
}
\examples{
## The data is one realisation from the the underlying DAG:
ex1.true.dag <- matrix(data=c(
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,
  1,1,0,0,0,0,0,0,0,0,
  1,0,1,1,0,0,0,0,0,0,
  0,1,1,1,0,0,0,0,0,0,
  0,0,1,0,1,0,0,0,0,0,
  0,0,1,0,0,0,1,0,0,0,
  0,0,1,1,0,0,0,0,0,0), ncol=10, byrow=TRUE)

colnames(ex1.true.dag) <- rownames(ex1.true.dag) <-
  c("b1","p1","g1","b2","p2","b3","g2","b4","b5","g3")
}
\keyword{datasets}
\keyword{internal}
