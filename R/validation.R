#' Make DAG of class "abnDag"
#'
#' @param dag matrix or formula specifying DAG.
#' @param data.df data frame
#' @param data.dists list of variable distributions.
#' @param ... additional arguments.
#'
#' @return abnDag object as list of dag, data.df, data.dists.
#'
#' Create a legitimate DAGs
#'
#' Create a legitimate DAG in the abn format.
#'
#' @usage
#' createAbnDag( dag, data.df = NULL, data.dists = NULL, ...)
#'
#' @param dag a matrix or a formula specifying a DAG, see \sQuote{Details}.
#' @param data.df named dataframe.
#' @param data.dists named list giving the distribution for each node in the network. If not provided it will be sample and returned.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' An object of class \code{class(abnDag)} contains a named matrix describing the DAG and possibly additional objects such as the associated distributions of the nodes.
#' If the dag is specified with a formula, either \code{data.df} or \code{data.dists} is required with the \code{.} quantifier.
#' If the dag is specified with an unnamed matrix and both \code{data.df} and \code{data.dists} are missing,  lower-case letters of the Roman alphabet are used for the node names.
#'
#' @return An object of class \code{abnDag} containing a named matrix and a named list giving the distribution for each node.
#' @export
#' @examples
#' dagFromFormula <- createAbnDag(dag = ~a+b|a,
#'                               data.df = data.frame("a"=1, "b"=1),
#'                               data.dists = list(a="binomial", b="gaussian"))
#' dagFromMatrix <- createAbnDag(dag = matrix(c(0,1,0,0), 2, 2),
#'                               data.df = data.frame("a"=1, "b"=1),
#'                               data.dists = list(a="binomial", b="gaussian"))
#' plot(dagFromMatrix)
#'
#' @keywords utilities internal
#' @concept DAG
createAbnDag <- function(dag,
                         data.df = NULL,
                         data.dists = NULL,
                         ...) {
  if (!is.null(data.dists)) {
    data.dists <- validate_dists(data.dists, returnDists = TRUE)
  }

  if (is.null(data.df)) {
    dag <- validate_abnDag(dag, data.df = data.dists, returnDag = TRUE) # Little hack to get the data distribution names in the correct format
  } else{
    if(is.data.frame(data.df)) {
      dag <- validate_abnDag(dag, data.df = data.df, returnDag = TRUE)
    } else {
      stop("'data.df' must be a data.frame")
    }
  }

  if (is.null(dimnames(dag))) {
    dag <- provideDimnames(dag, base = list(letters))
    validate_abnDag(dag)
  }

  out <- list(dag = dag,
              data.df = data.df,
              data.dists = data.dists)
  class(out) <- "abnDag"

  return(out)
}

#' Check for valid distribution
#'
#' The distribution names must match `inla() family=''`.
#' Similar to `get.var.types()`, mainly different in output.
#'
#' @param data.dists list of variable distributions.
#' @param returnDists if TRUE (default) returns the same list as provided.
#' @param ... additional arguments.
#'
#' @return either TRUE/FALSE or list of variable distributions as provided as input.
#' @keywords internal
validate_dists <- function(data.dists, returnDists=TRUE,...) {

  name <- names(data.dists)
  if (is.null(name)) stop("Node distribution has to be a named object.")
  if( is.list( data.dists))       data.dists <- unlist( data.dists)

  choices <- c("poisson","binomial","gaussian","multinomial")
  data.dists <- choices[pmatch(tolower(data.dists ), choices, duplicates.ok=TRUE)]
  if (any(is.na(data.dists ))){
    ## each variable must have one of the allowed.dists
    e <- simpleError(paste("Each variable must have a distribution of either", paste(choices, collapse = ", ")))
    stop(e)
  }
  names(data.dists ) <- name


  if( returnDists) return( as.list( data.dists)) else return( TRUE)

}

#' Check for valid DAG of class `abnDag`
#'
#' Beside some basic checks, this function also checks for square matrix with no undirected cycles (trivial cycles) and
#' for no undirected cycles.
#'
#' Similar to `check.valid.dag()`.
#'
#' @param dag dag is either a formula, a matrix  or an object of class 'abnDag'
#' @param data.df data frame
#' @param returnDag if TRUE (default) returns DAG as matrix.
#' @param ... additional arguments.
#'
#' @return Either TRUE/FALSE or DAG as matrix.
#' @keywords internal
validate_abnDag <- function( dag, data.df=NULL, returnDag=TRUE, ...) {
    # we already have a valid container. can be used to extract...
  if (inherits(x = dag, what = "abnDag"))  dag <- dag$dag


  # case of formula
  if  (inherits(x = dag, what = "formula")) {
    if (is.null( data.df))
      stop( 'DAG specification with formula requires a named data frame or named vector')

    name <- if ( is.matrix( data.df)) colnames( data.df) else names( data.df)
    if (is.null( name))
      stop( 'Improperly named object "data.df" for DAG specification')

    dag <- formula_abn(f = dag, name = name)
  }   # proceed checking!!

  # case of matrix

  if ( is.matrix( dag)) {
    dimm <- dim( dag)
    if (dimm[1] != dimm[2])   stop("DAG matrix is not square")
    if (any(diag(dag)!=0))  stop("DAG matrix contains trivial cycles (nonzero values on diagonal)")

    if (!is.null(data.df))  {    # if data.df given we take over the names.
      name <- if ( is.matrix( data.df)) colnames( data.df) else names( data.df)

      if(length(name) != dimm[1])  stop("DAG matrix not coherent with names")
      colnames(dag) <- rownames(dag) <- name
    } else {
      if (any(colnames(dag)!=rownames(dag)))  stop("DAG matrix with incoherent row-/colnames")
    }

    res <- .Call("checkforcycles", as.integer(dag), dimm[1], PACKAGE = "abn")
    if (res!=0) stop("DAG contains at least one cycle.")

    if( returnDag) return( dag) else return( TRUE)
  }   else {
    stop( "DAG specification with should be via formula or matrix")
  }
}
