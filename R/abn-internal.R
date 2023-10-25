#' abn Version Information
#'
#' \code{abn.version()} provides detailed information about the running version of \pkg{abn}
#' or the \pkg{abn} components.
#'
#' @param what detailed information about the version of \pkg{abn} or the system (see returns).
#'
#' @return \code{abn.version(what = "system")} is a list with character-string components
#' \describe{
#' \item{R}{\code{R.version.string}}
#' \item{abn}{essentially \code{abn.version$version.string}}
#' \item{GSL, JAGS, INLA}{version numbers thereof}
#' }
#'
#' \code{abn.version(what = "abn")} is a list with character-string components
#' \describe{
#' \item{status}{the status of the version (e.g., \code{"beta"})}
#' \item{major}{the major version number}
#' \item{minor}{the minor version number}
#' \item{year}{the year the version was released}
#' \item{month}{the month the version was released}
#' \item{day}{the day the version was released}
#' \item{version.string}{a \code{character} string concatenating
#' the info above, useful for plotting, etc.}
#' }
#'
#' \code{abn.version} is a list of class \code{"simple.list"} which has a \code{print} method.
#' @seealso \code{\link[base]{R.version}}
#' @export
#'
#' @examples
#' abn.version()$version.string
#' \dontrun{
#'   abn.version("system")
#' }
abn.version <- function(what=c('abn','system')) {
  what <- match.arg(what)
  if (what %in% 'system') {

    list(R=R.version.string,
         abn=substr(abn.version()$version.string, 13, 32),
         gsl=ifelse(R.version$os=="linux-gnu", system('gsl-config --version', intern = TRUE), "NA (?)"),
         JAGS=rjags::jags.version(),
         INLA=ifelse(requireNamespace("INLA", quietly = TRUE),
           INLA::inla.version("version"), "not available")
    )

  } else {
    release <- utils::packageDescription("abn",field="Version")
    date <- utils::packageDescription("abn",field="Date")
    list(status="",
        major=sub("-","",substr(release,1,4)),
        minor=substr(sub("-","",substr(release,5,7)),1,1),
         year=substr(date,1,4),
         month=substr(sub("20..-","",date),1,2),
         day=sub("20..-..-","",date),
         version.string= paste("abn version ",
                             utils::packageDescription("abn",field="Version")," (",
                             utils::packageDescription("abn",field="Date"),")",sep="")
  )
  }
}

#' @title Prints start up message
#' @examples
#' library(abn)
#' @keywords internal
".onAttach" <- function (lib, pkg) {
  packageStartupMessage(abn.version()$version.string," is loaded.\nTo cite the package 'abn' in publications call: citation('abn').")
}

#' @title Recursive string splitting
#' @description Internal function that call multiple times strsplit() and remove space
#' @keywords internal
#' @export
strsplits <- function(x, splits, ...) {
    for (split in splits) {
        x <- unlist(strsplit(x, split, ...))
    }
    x <- gsub(" ", "", x, fixed = TRUE)  #remove space
    return(x[!x == ""])  # Remove empty values
}

#' @title Formula to adjacency matrix
#' @description Internal function that produce a square matrix length(name) with \eqn{0,1} depending on f.
#' f have to start with ~ terms are entries of name terms are separated by + term1 | term2 indicates
#' col(term1) row(term2) puts a 1 term1 | term2:term3: ... : is used as a sep . = all terms in name
#' @keywords internal
#' @export
formula_abn <- function(f, name) {

    name_orignial <- name

    f <- as.character(f)

    ## tests for consistence ---------------------------------------------------------------------- start as a formula
    if (!grepl("~", f[1], fixed = T)) {
        stop("DAG specifications should start with a ~")
    }

    ## transformation name + or | or : or . or name to name_name
    if (sum((c("+", "|", ":", ".") %in% unlist(strsplit(name, split = c(""))))) != 0) {
        for (i in 1:length(name)) {
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c("+")) != 0) {
                f[[2]] <- gsub(name[i], gsub("+", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub("+", "_", name[i], fixed = TRUE)
            }
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c("|")) != 0) {
                f[[2]] <- gsub(name[i], gsub("|", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub("|", "_", name[i], fixed = TRUE)
            }
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c(":")) != 0) {
                f[[2]] <- gsub(name[i], gsub(":", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub(":", "_", name[i], fixed = TRUE)
            }
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c(".")) != 0) {
                f[[2]] <- gsub(name[i], gsub(".", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub(".", "_", name[i], fixed = TRUE)
            }
        }
    }

    ## collapse name
    name.c <- paste(name, collapse = ":")
    ## Split by terms
    f.p <- strsplit(x = f[[2]], split = "+", fixed = TRUE)

    ## nothing more than name variable in the dag formula
    tmp.test <- strsplits(x = f[[2]], splits = c("+", "|", ":", "."), fixed = TRUE)
    if (sum(!(tmp.test %in% name)) != 0) {
        stop("DAG formulation contains some variables not in provided names")
    }
    ## End of tests for consistence ----------------------------------------------------------------

    ## creat the void matrix
    out <- matrix(data = 0, nrow = length(name), ncol = length(name))

    ## delete all spaces
    f.p <- gsub(" ", "", f.p[[1]], fixed = TRUE)

    ## replace '.' by all names
    f.p.completed <- gsub(".", name.c, f.p, fixed = TRUE)

    ## atomization of left term


    ## contruction of the output matrix
    for (i in 1:length(f.p)) {
        tmp <- f.p.completed[i]

        ## forget unique terms -> test for |
        if (grepl("|", tmp, fixed = TRUE)) {

            ## split wrt |
            tmp.p <- strsplit(x = tmp, split = "|", fixed = TRUE)

            ## test for multiple terms and contruction of the list first term
            if (grepl(":", tmp.p[[1]][1])) {
                tmp.p.p.1 <- strsplit(x = tmp.p[[1]][1], split = ":", fixed = TRUE)
            }
            if (!grepl(":", tmp.p[[1]][1])) {
                tmp.p.p.1 <- tmp.p[[1]][1]
            }

            ## test for multiple terms and contruction of the list second term
            if (grepl(":", tmp.p[[1]][2])) {
                tmp.p.p.2 <- strsplit(x = tmp.p[[1]][2], split = ":", fixed = TRUE)
            }
            if (!grepl(":", tmp.p[[1]][2])) {
                tmp.p.p.2 <- tmp.p[[1]][2]
            }

            ## loop over the
            for (j in 1:length(tmp.p.p.1[[1]])) {
                for (k in 1:length(tmp.p.p.2[[1]])) {
                  ## update of matrix
                  out[grep(tmp.p.p.1[[1]][j], name), grep(tmp.p.p.2[[1]][k], name)] <- 1

                }
            }
        }

    }

    ## avoid auto dependance
    diag(out) <- 0

    ## only 0 and 1
    out[out > 1] <- 1

    ## naming
    colnames(out) <- name_orignial
    rownames(out) <- name_orignial
    ## output
    return(out)
}

#' Set of simple commonsense validity checks on the data.df and data.dists arguments
#'
#' @param data.df data frame with names corresponding to variable/node names.
#' @param data.dists list specifying each columns distribution type. Names correspond to column names and values must be one of "gaussian", "binomial", "poisson", "multinomial".
#' @param group.var not yet implemented
#'
#' @return list of sums of each distribution types (abbreviated) as names.
#' @importFrom stats complete.cases
check.valid.data <- function(data.df = NULL, data.dists = NULL, group.var = NULL) {

    ## check data is in a data.frame
    if (!is.data.frame(data.df)) {
        stop("The data must be in a data.frame")
    }

    ## check data for missing values
    if (sum(complete.cases(data.df)) != dim(data.df)[1]) {
        stop("The data contains missing values! These must be removed.")
    }

    ## check that distributions are in a list
    if (!is.list(data.dists)) {
        stop("data.dist must be a list")
    }

    if (!is.null(group.var)) {
        ## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY
      data.df <- check.valid.groups(group.var = group.var, data.df = data.df)[["data.df"]]
    }

    if (length(names(data.dists)) != length(names(data.df))) {
        stop("data.dists must have named entries")
    }

    ## check names in list are in correct order and exact match to data.frame
    for (i in 1:dim(data.df)[2]) {
        if (names(data.dists)[i] != names(data.df)[i]) {
            stop("names in list must match names in data.frame")
        }
    }

    ## check names of likelihood function are valid
    validate_dists(data.dists = data.dists, returnDists = FALSE)

    binomial.vars.indexes <- NULL
    poisson.vars.indexes <- NULL
    gaussian.vars.indexes <- NULL
    multinomial.vars.indexes <- NULL

    ## check that data is consistent with distribution given for each variable
    for (i in 1:dim(data.df)[2]) {
        cur.var <- data.df[, i]
        if (data.dists[[i]] == "gaussian") {
            if (is.factor(cur.var)) {
                cat((names(data.df)[i]), "is invalid - it must not be a factor.\n")
                stop("")
            }
            if (length(unique(cur.var)) <= 2) {
                cat((names(data.df)[i]), "is invalid as it has two or less unique values!\n")
                stop("")
            }
            gaussian.vars.indexes <- c(gaussian.vars.indexes, i)
        }
        if (data.dists[[i]] == "binomial") {
            if (!is.factor(cur.var)) {
                cat((names(data.df)[i]), "is invalid - it must be a factor\n")
                stop("")
            }
            if (length(unique(cur.var)) < 2) {
                cat((names(data.df)[i]), "is invalid as it must be binary with both cases being observed.\n")
                stop("")
            }
            if (length(unique(cur.var)) > 2) {
                cat((names(data.df)[i]), "is invalid as it must be binary. Multi-category variables should be split into separate binary variables or defined as multinomial distributed.\n")
                stop("")
            }
            binomial.vars.indexes <- c(binomial.vars.indexes, i)
        }

        if (data.dists[[i]] == "poisson") {
            if (is.factor(cur.var)) {
                cat((names(data.df)[i]), "is invalid - it must not be a factor\n")
                stop("")
            }
            if (length(unique(cur.var)) <= 2) {
                cat((names(data.df)[i]), "is invalid as it has two or less unique values!\n")
                stop("")
            }
            poisson.vars.indexes <- c(poisson.vars.indexes, i)
        }

        if (data.dists[[i]] == "multinomial") {
          if (!is.factor(cur.var)) {
            cat((names(data.df)[i]), "is invalid - it must be a factor\n")
            stop("")
          }
          if (length(unique(cur.var)) <= 2) {
            cat((names(data.df)[i]), "is invalid as it has two or less unique values! Consider it as binary variable.\n ")
            stop("")
          }
          multinomial.vars.indexes <- c(multinomial.vars.indexes, i)
        }
    }
    ## return the indexes of any binary variables
    return(list(gaus = gaussian.vars.indexes, bin = binomial.vars.indexes, pois = poisson.vars.indexes, mult = multinomial.vars.indexes))

}  #end of check.valid.data()

#' Set of simple commonsense validity checks on the directed acyclic graph definition matrix
#'
#' @param dag Named square matrix or a formula statement specifying a directed acyclic graph.
#' If NULL an empty network is assumed.
#' @param data.df data frame with names corresponding to variable/node names.
#' @param is.ban.matrix Diagonals and cycles are not checked for ban-matrices. Defaults to FALSE.
#' @param group.var not yet implemented
#'
#' @return dag as named square matrix
check.valid.dag <- function(dag = NULL, data.df = NULL, is.ban.matrix = FALSE, group.var = NULL) {
  if (!is.null(data.df) && !inherits(data.df, "data.frame")){
    stop("Invalid argument for data.df provided. Must be NULL or of class data.frame.")
  }

    if (!is.null(group.var)) {
        ## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY

        if (is.null(data.df)) stop("When specifying 'group.var', 'data.df' argument is required as well.")
        data.df <- data.df[, -which(names(data.df) == group.var)]
    }


    ## if dag null then create unlimited - empty - network want ban matrix
    if (is.null(dag)) {
        dag <- matrix(rep(0, dim(data.df)[2]^2), ncol = dim(data.df)[2])
        ## names must be set
        colnames(dag) <- rownames(dag) <- names(data.df)
        return(dag)
    }

    if (!is.matrix(dag)) {
        if (is.null(data.df)){
          stop("If dag is provided as formula, data.df must be specified.")
        } else if(grepl("~", as.character(dag)[1], fixed = T)) {
          # if provided as formula, convert to matrix representation and do all the tests (e.g. acyclicity!)
            dag <- formula_abn(f = dag, name = names(data.df))
        } else {
            stop("'dag' specification must either be a matrix or a formula expression.")
        }
    } else {
        if (dim(dag)[1] != dim(dag)[2]){
          stop("Matrix 'dag' is not square.")
        }
    }

    ## check data for missing names
    if (is.null(colnames(dag)) || is.null(rownames(dag))) {
        if (!is.null(data.df)) {
            if (dim(dag)[1] != dim(data.df)[2]) {
              stop("'dag' as dimension inconsistent with columns of 'data.df'")
            }
          colnames(dag) <- rownames(dag) <- names(data.df)
        } else {
            stop("'dag' must be a matrix with both row and column names set or a named dataset has to be provided")
        }
    }
    ## check dimension
    if (!is.null(data.df)) {
      if (dim(dag)[1] != dim(data.df)[2] || dim(dag)[2] != dim(data.df)[2]) {
         stop("'dag' as dimension inconsistent with data.df")
      }
    }

  ## check equal order of row and column names in dag
  if(any(colnames(dag) != rownames(dag))){
    stop("dag must be a symmetric, named matrix: Row and column names must be in the same order.")
  }

    ## check binary
    for (i in 1:dim(dag)[1]) {
       for (j in 1:dim(dag)[2]) {
          if ( abs(dag[i, j]) > 1e-8 && abs(dag[i, j]-1) > 1e-8)      stop("'dag' must comprise only 1's or 0's")
       }
    }
    # if (any( c( (dag != 0) && (dag != 1))))                stop("'dag' must comprise only 1's or 0's")


    ## check diagnonal and cycles - but ignore these checks for a ban matrices
    if (!is.ban.matrix) {
      if(any(diag(dag) != 0)){
        stop("'dag' is not a valid DAG - a child cannot be its own parent!")
      }
      ## coerce to int for sending to C$ number of cols (or rows)
      ## this creates one long vector - filled by cols from dag = same as internal C reprentation so fine.
      res <- .Call("checkforcycles", as.integer(dag),  dim(dag)[1], PACKAGE = "abn")
      if (res!=0){
        stop("'dag' contains at least one cycle.")
      }
    }
  return( dag)
}


#' Set of simple checks on the given parent limits
#'
#' @param data.df data frame
#' @param max.parents single integer for one overall max parent limit.
#' A list with names corresponding to variable/column names of `data.df` and individual parent limits.
#' NULL for no parent limit restriction(s).
#' @param group.var not yet implemented
#'
#' @return numeric vector of max number of parents per variable
check.valid.parents <- function(data.df = NULL, max.parents = NULL, group.var = NULL) {
  ## Stop if data.df is not provided
  if (is.null(data.df)){
    stop("`data.df` is not provided.")
  }

  ## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY
  if (!is.null(group.var)) {
    data.df <- check.valid.groups(group.var = group.var, data.df = data.df)[["data.df"]]
  } else {
    # TODO: Take care of group.var
  }

  if (is.numeric(max.parents)){
    # if a constant then make integer vector
    if(length(max.parents) == 1) {
      # check if no value is larger than the number of variables in data.df
      if (max.parents > dim(data.df)[2]){
        stop("`max.parents` is larger than the total number of variables.")
      } else {
        return(as.integer(rep(max.parents, dim(data.df)[2])))
      }
    } else if (length(max.parents > 1)){
      # check if no value is larger than the number of variables in data.df
      if (any(max.parents > dim(data.df)[2])){
        stop("`max.parents` has values that are larger than the total number of variables.")
      } else {
          return(max.parents)
      }
    } else {
      stop("max.parents must be either a named list or an integer vector with their size corresponding to the number of variables. Alternatively, it can be a single integer.")
    }
  } else if (is.list(max.parents)){
    ## if a list must be named list with names as in original data
    if(any(mapply(is.null, max.parents)) || any(mapply(is.na, max.parents))){
      stop("Values of max.parents list should not be empty (NULL or NA).")
    } else if (length(max.parents) != dim(data.df)[2]){
      stop("Length of max.parents (", length(max.parents), ") is not equal to the number of variables. Provide for each variable an individual number of max.parents.")
    } else if (length(max.parents) == dim(data.df)[2]) {
      if (any(!mapply(is.numeric, max.parents))) {
        stop("max.parents is not valid - must be numeric")
      }

      for (i in 1:dim(data.df)[2]) {
        if (names(max.parents)[i] != names(data.df)[i]) {
          stop("names in max.parents list must match names in data.frame data.df")
        } else if (max.parents[[i]] > dim(data.df)[2]){
          stop(paste("`max.parents` has values that are larger (", max.parents[i], ") than the total number of variables (",  dim(data.df)[2], ")."))
        }
      }

      max.parents.int <- unlist(max.parents, use.names = FALSE)
      if (length(max.parents.int) != dim(data.df)[2]) {
        stop("max.parents list has wrong length")
      }

      max.parents.int <- as.integer(max.parents.int)
      return(max.parents.int)
    }
  } else if (is.null(max.parents)) {
    ## if NULL, return integer vector of max possible values
    return(as.integer(rep(max(dim(data.df)[2]), dim(data.df)[2])))
  } else {
    stop("max.parents must be either numeric, a list or NULL.")
  }

  # Raise error if situation was not catched above.
  stop("'max.parents' is not valid: length data: ",dim(data.df)[2],
       ", length max.parents: ",length(max.parents))
}

#' Set of simple checks on the list given as parent limits
#'
#' @param data.df data frame
#' @param which.nodes a vector giving the column indices of the variables to be included, if ignored all variables are included.
#' @param group.var not yet implemented
#'
#' @return integer vector of column indexes of valid nodes in data.df
check.which.valid.nodes <- function(data.df = NULL, which.nodes = NULL, group.var = NULL) {
    ## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY
    if (!is.null(group.var)) {
      data.df <- check.valid.groups(group.var = group.var, data.df = data.df)[["data.df"]]
    }
    ## if null then assume ALL nodes
    if (is.null(which.nodes)) {
        which.nodes <- 1:dim(data.df)[2]
        return(as.integer(which.nodes))
    }

    if (is.numeric(which.nodes) && max(which.nodes) <= dim(data.df)[2] && min(which.nodes) >= 1 && length(which.nodes) <= dim(data.df)[2]) {
        return(as.integer(which.nodes))
    } else {
        stop("which.nodes is invalid")
    }

}

#' Simple check on the grouping variable
#'
#' @param group.var Name of grouping variable.
#' @param data.df data frame of all variables including the variable specified in `group.var` as factor.
#' @param cor.vars Name(s) of variables to which the grouping should be applied to.
#' @param verbose when TRUE additional information is printed. Defaults to FALSE.
#'
#' @return list with data.df, indexes of variables to which the grouping should be applied to and the associated group for each observation as integer.
check.valid.groups <- function(group.var=NULL, data.df=NULL, cor.vars=NULL, verbose = FALSE) {
  # No data no checks.
  if (is.null(data.df)){
    stop("No data.df provided.")
  }

  ## have no cor.vars, take all but group.var.
  if (is.null(cor.vars)){
    if (verbose){
      warning("No cor.vars specified. Using all but group.var instead.")
    }
    cor.vars <- names(data.df[, which(names(data.df) != group.var)])
  } else if (is.null(group.var)){
    stop("If cor.vars is given, group.var must be specified too!")
  }

  ## have no groups so just return dummy values
  if (is.null(group.var)){
    return(list(data.df = data.df, grouped.vars = as.integer(c(-1)), group.ids = as.integer(rep(-1, dim(data.df)[1]))))
  }

  ## Check group.var
  if (!(is.character(group.var) && (length(group.var) == 1))) {
    stop("name of group variable is not a character?!")
  }
  if (!length(which(group.var %in% names(data.df) == TRUE))) {
    stop("name of group variable does not match any of those in data.df")
  }

  ## Check cor.var
  if (!(is.character(cor.vars))) {
    stop("name of cor.var is not a character?!")
  }
  if (!length(which(cor.vars %in% names(data.df) == TRUE))) {
    stop("name of cor.vars does not match any of those in data.df")
  }
  if (group.var %in% cor.vars) {
    stop("group.var is among the cor.vars.") # TODO: consider to relax this.
  }

  ## get group id data
  group.var.vals <- data.df[, group.var]
  ## drop the group variable from original data.frame and overwrite
  data.df <- data.df[, -which(names(data.df) == group.var)]

  ## have groups so some checks
  if (is.factor(group.var.vals) && length(group.var.vals) == dim(data.df)[1] && length(unique(group.var.vals)) > 1) {
    ## is factor and of correct length and at least two groups
  } else {
    stop("grouping variable must be: i) a factor; ii) same number of observations as in data.df; and iii) contain at least two different observations (levels)")
  }

  ## get group memberships in terms of ints
  group.var <- as.integer(group.var.vals)

  ## now find out which variables the grouping is to be applied to
  var.noms <- names(data.df)
  if (length(which(cor.vars %in% var.noms == TRUE)) != length(cor.vars)) {
    stop("variables in cor.vars do not match those in data.df")
  }

  if (max(table(cor.vars)) > 1) {
    stop("have repeated variables in cor.vars!")
  }

  ## to get to here group.var must be ok and also variable names so return integer code for the variables get the index in names(data.df) for each variable and then sort into order
  cor.var.indexes <- as.integer(sort(match(cor.vars, var.noms)))

  return(list(data.df = data.df, grouped.vars = cor.var.indexes, group.ids = group.var))
}

#' Simple check on the control parameters
#'
#' @param control list of control arguments with new parameters supplied to \code{\link{buildScoreCache}} or \code{\link{fitAbn}}.
#' @param method "bayes" or "mle" strategy from argument \code{method=...} in \code{\link{buildScoreCache}} or \code{\link{fitAbn}}. Defaults to "bayes".
#' @param verbose when TRUE additional information is printed. Defaults to FALSE.
#'
#' @return list with all control arguments with respect to the method but with new values.
check.valid.buildControls <- function(control, method = "bayes", verbose = FALSE) {
  ctrl.basic <- build.control(method = method)
  if (is.null(control)) {
    return(build.control(method = method))
  } else if (!is.null(control)){
    ctrl.new <- control
  } else {
    stop("Invalid 'control' argument.")
  }

  # check type of control
  if(!is.list(ctrl.new)) {
    stop("Control arguments must be provided as named list.")
  }

  # check if keys are ok
  allowed_list_names <- c("method",
                          "max.mode.error",
                          "mean",
                          "prec",
                          "loggam.shape",
                          "loggam.inv.scale",
                          "max.iters",
                          "epsabs",
                          "error.verbose",
                          "trace",
                          "epsabs.inner",
                          "max.iters.inner",
                          "finite.step.size",
                          "hessian.params",
                          "max.iters.hessian",
                          "max.hessian.error",
                          "factor.brent",
                          "maxiters.hessian.brent",
                          "num.intervals.brent",
                          "n.grid",
                          "ncores",
                          "max.irls",
                          "tol",
                          "tolPwrss",
                          "check.rankX",
                          "check.scaleX",
                          "check.conv.grad",
                          "check.conv.singular",
                          "check.conv.hess",
                          "xtol_abs",
                          "ftol_abs",
                          "trace.mblogit",
                          "catcov.mblogit",
                          "epsilon",
                          "seed")
  if(any(!(names(ctrl.new) %in% allowed_list_names))) {
    stop("Unknown control parameter(s).")
  } else if(any(!(names(ctrl.new) %in% names(build.control(method=method))))) {
    warning(paste("Control parameters provided that are not used with method", method, "are ignored."))
    # ctrl <- ctrl.basic[which(!(names(build.control(method=method)) %in% names(ctrl.new)))] # ignore them further down in collecting list for return
  }

  # TODO: Add more checks here for the individual control parameters.
  # check catcov.mblogit
  possible_catcov.mblogit <- c("free", "diagonal", "single")
  if (!is.null(ctrl.new[["catcov.mblogit"]])) {
    if (!(ctrl.new[["catcov.mblogit"]] %in% possible_catcov.mblogit)) {
      stop(paste("'catcov.mblogit' must be one of", deparse(possible_catcov.mblogit)))
    }
  }
  # check max.mode.error
  if (!is.null(ctrl.new[["max.mode.error"]])) {
    if (!((ctrl.new[["max.mode.error"]] >= 0) && (ctrl.new[["max.mode.error"]] <= 100))) {
      stop("'max.mode.error' is a % and must be [0,100]!")
    }
  }
  # check ncores
  if (!is.null(ctrl.new[["ncores"]])) {
    if (method != "mle") {
      warning("Multithreading is currently only implemented for method='mle'. I'm ignoring 'ncores' and continue with a single core.")
      ctrl.new[["ncores"]] <- 1
    } else {
      # Prepare multithreading
      if (ctrl.new[["ncores"]] == -1) {
        # all but one
        ctrl.new[["ncores"]] <-  parallel::detectCores() - 1   # if ncores==0 (here or set), single threaded.
        if(verbose){message("Running in parallel with ", ctrl.new[["ncores"]], " cores.")}
      } else if (ctrl.new[["ncores"]] > 1) {
        ctrl.new[["ncores"]] <- min(ctrl.new[["ncores"]], parallel::detectCores())  # restrict in case of overoptimisitic setting.
        if(verbose){message("Running in parallel with ", ctrl.new[["ncores"]], " cores.")}
      } else if (ctrl.new[["ncores"]] == 1 | ctrl.new[["ncores"]] == 0) {
        ctrl.new[["ncores"]] <- 1
        if(verbose){message("Running in single core mode.")}
      } else {
        stop(paste("Argument 'ncores' from build.control(ncores=...) is invalid. It must be larger or equal -1 and smaller or equal", parallel::detectCores()))
      }
    }
  }
  # check seed
  if (!is.null(ctrl.new[["seed"]])) {
    if ((!inherits(ctrl.new[["seed"]], "integer") | ctrl.new[["seed"]] < 0)) {
      stop("'seed' must be a non-negative integer.")
    }
  }
  # Collect control list for return
  ctrl <- ctrl.basic
  for (i in names(ctrl.new)){
    if (!is.null(ctrl[[i]])) {
      # if new parameter is used by the specific method, update parameter
      ctrl[[i]] <- ctrl.new[[i]]
    } else {
      # ignore parameters that are not used by current method
    }
  }
  return(ctrl)
}

#' @inherit check.valid.buildControls
check.valid.fitControls <- function(control, method = "bayes", verbose = FALSE) {
  ctrl.basic <- fit.control(method = method)
  if (is.null(control)) {
    return(fit.control(method = method))
  } else if (!is.null(control)){
    ctrl.new <- control
  } else {
    stop("Invalid 'control' argument.")
  }

  # check type of control
  if(!is.list(ctrl.new)) {
    stop("Control arguments must be provided as named list.")
  }

  # check if keys are ok
  allowed_list_names <- c("method",
                          "max.mode.error",
                          "mean",
                          "prec",
                          "loggam.shape",
                          "loggam.inv.scale",
                          "max.iters",
                          "epsabs",
                          "error.verbose",
                          "trace",
                          "epsabs.inner",
                          "max.iters.inner",
                          "finite.step.size",
                          "hessian.params",
                          "max.iters.hessian",
                          "max.hessian.error",
                          "factor.brent",
                          "maxiters.hessian.brent",
                          "num.intervals.brent",
                          "min.pdf",
                          "n.grid",
                          "std.area",
                          "marginal.quantiles",
                          "max.grid.iter",
                          "marginal.node",
                          "marginal.param",
                          "variate.vec",
                          "ncores",
                          "max.irls",
                          "tol",
                          "tolPwrss",
                          "check.rankX",
                          "check.scaleX",
                          "check.conv.grad",
                          "check.conv.singular",
                          "check.conv.hess",
                          "xtol_abs",
                          "ftol_abs",
                          "trace.mblogit",
                          "catcov.mblogit",
                          "epsilon",
                          "seed")
  if(any(!(names(ctrl.new) %in% allowed_list_names))) {
    stop("Unknown control parameter(s).")
  } else if(any(!(names(ctrl.new) %in% names(fit.control(method=method))))) {
    warning(paste("Control parameters provided that are not used with method", method, "are ignored."))
  }

  # TODO: Add more checks here for the individual control parameters.
  # check catcov.mblogit
  possible_catcov.mblogit <- c("free", "diagonal", "single")
  if (!is.null(ctrl.new[["catcov.mblogit"]])) {
    if (!(ctrl.new[["catcov.mblogit"]] %in% possible_catcov.mblogit)) {
      stop(paste("'catcov.mblogit' must be one of", deparse(possible_catcov.mblogit)))
    }
  }
  # Check max.grid.iter
  if(!is.null(ctrl.new[["max.grid.iter"]])) {
    if(is.na(ctrl.new[["max.grid.iter"]])) {
      stop("'max.grid.iter' cannot be NA!")
    } else {
      if(!is.null(ctrl.new[["variate.vec"]])) {
        ctrl.new[["max.mode.error"]] <- 0 # if user supplied grid then must use C
        ctrl.new[["std.area"]] <- FALSE
        ctrl.new[["n.grid"]] <- NULL;
      }
    }
  }

  # check max.mode.error
  if (!is.null(ctrl.new[["max.mode.error"]])) {
    if (!((ctrl.new[["max.mode.error"]] >= 0) && (ctrl.new[["max.mode.error"]] <= 100))) {
      stop("'max.mode.error' is a % and must be [0,100]!")
    }
  }
  # check ncores
  if (!is.null(ctrl.new[["ncores"]])) {
    if (method != "mle") {
      warning("Multithreading is currently only implemented for method='mle'. I'm ignoring 'ncores' and continue with a single core.")
      ctrl.new[["ncores"]] <- 1
    } else {
      # Prepare multithreading
      if (ctrl.new[["ncores"]] == -1) {
        # all but one
        ctrl.new[["ncores"]] <-  parallel::detectCores() - 1   # if ncores==0 (here or set), single threaded.
        if(verbose){message("Running in parallel with ", ctrl.new[["ncores"]], " cores.")}
      } else if (ctrl.new[["ncores"]] > 1) {
        ctrl.new[["ncores"]] <- min(ctrl.new[["ncores"]], parallel::detectCores())  # restrict in case of overoptimisitic setting.
        if(verbose){message("Running in parallel with ", ctrl.new[["ncores"]], " cores.")}
      } else if (ctrl.new[["ncores"]] == 1 | ctrl.new[["ncores"]] == 0) {
        ctrl.new[["ncores"]] <- 1
        if(verbose){message("Running in single core mode.")}
      } else {
        stop(paste("Argument 'ncores' from fit.control(ncores=...) is invalid. It must be larger or equal -1 and smaller or equal", parallel::detectCores()))
      }
    }
  }
  # check seed
  if (!is.null(ctrl.new[["seed"]])) {
    if ((!inherits(ctrl.new[["seed"]], "integer") | ctrl.new[["seed"]] < 0)) {
      stop("'seed' must be a non-negative integer.")
    }
  }
  # Collect control list for return
  ctrl <- ctrl.basic
  for (i in names(ctrl.new)){
    if (!is.null(ctrl[[i]])) {
      # if new parameter is used by the specific method, update parameter
      ctrl[[i]] <- ctrl.new[[i]]
    } else {
      # ignore parameters that are not used by current method
    }
  }
  return(ctrl)
}


#' Create ordered vector with integers denoting the distribution
#'
#' gaussian = 1, binomial = 2, poisson = 3, multinomial = 4
#'
#' @param data.dists list specifying each columns distribution type. Names correspond to column names and values must be one of "gaussian", "binomial", "poisson", "multinomial".
#'
#' @return numeric encoding of distribution corresponding to its list element number in `data.dists`.
get.var.types <- function(data.dists = NULL) {
  store <- rep(NA, length(data.dists))

  for (i in 1:length(data.dists)) {
    if (data.dists[[i]] == "binomial") {
      store[i] <- 1
    }
    if (data.dists[[i]] == "gaussian") {
      store[i] <- 2
    }
    if (data.dists[[i]] == "poisson") {
      store[i] <- 3
    }
    if (data.dists[[i]] == "multinomial") {
      store[i] <- 4
    }
  }

  # Check if all distributions could be encoded
  if (sum(is.na(store)) > 0){
    # List elements with failed distribution encoding
    failed_idx <- which(is.na(store))
    stop("Unknown distribution type(s): \n", as.data.frame(unlist(data.dists[failed_idx]))) # a bit unfortunate print out formatting.
  } else {
    return(store)
  }
}

#' @title tidy up cache
#' @keywords internal
#' @export
tidy.cache <- function(thecache) {
    if (!is.null(thecache[["error.indexes"]])) {
        error.combs <- thecache[["error.indexes"]]
        corrected <- list()
        corrected[["children"]] <- thecache[["children"]][-error.combs]
        corrected[["node.defn"]] <- thecache[["node.defn"]][-error.combs, ]
        corrected[["mlik"]] <- thecache[["mlik"]][-error.combs]
        ## return new cache with appropriate nodes removed
        return(corrected)
    }

}

