% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getmarginals.R
\name{find.next.left.x}
\alias{find.next.left.x}
\alias{find.next.right.x}
\title{Find next X evaluation Point}
\usage{
find.next.left.x(mat.xy, g.max, g.factor, x.delta, max.fact.delta)

find.next.right.x(mat.xy, g.max, g.factor, x.delta, max.fact.delta)
}
\arguments{
\item{mat.xy}{matrix}

\item{g.max}{integer. See Details.}

\item{g.factor}{integer. See Details.}

\item{x.delta}{integer. See Details.}

\item{max.fact.delta}{integer. See Details.}
}
\value{
integer

integer
}
\description{
Attempt to find the next x evaluation point using spline extrapolation traversing left from mode.
}
\details{
if new x point is more than a factor max.fact.delta (e.g. 0.2) from last evaluated point then stop here
\code{cat("evaluating node ",nodeid,": parameter:",paramid," at betafixed=",betafixed," with gvalue=",gvalue,"\n",sep="");}
find the next x value left which differs from the max. gvalue by at least a factor of g.factor, searching in step sizes of
x.delta subject to the constraint that if we move more than max.fact.delta*last.x then we evaluate here. Avoids big steps.
}
\section{Functions}{
\itemize{
\item \code{find.next.right.x()}: Attempt to find the next x evaluation point using spline extrapolation traversing right from mode.

}}
\keyword{internal}
