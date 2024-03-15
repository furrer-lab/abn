#' Fit a given regression using INLA
#'
#' Internal wrapper to INLA and are called from \code{fitAbn.bayes} and \code{buildScoreCache.bayes}.
#'
#' @param child.loc index of current child node.
#' @param dag.m.loc dag as matrix.
#' @param data.df.loc data df,
#' @param data.dists.loc list of distributions.
#' @param ntrials.loc \code{rep(1,dim(data.df)[1])}.
#' @param exposure.loc \code{rep(1,dim(data.df)[1])}.
#' @param compute.fixed.loc TRUE.
#' @param mean.intercept.loc the prior mean for all the Gaussian additive terms for each node. INLA argument \code{control.fixed=list(mean.intercept=...)} and \code{control.fixed=list(mean=...)}.
#' @param prec.intercept.loc the prior precision for all the Gaussian additive term for each node. INLA argument \code{control.fixed=list(prec.intercept=...)} and \code{control.fixed=list(prec=...)}.
#' @param mean.loc the prior mean for all the Gaussian additive terms for each node. INLA argument \code{control.fixed=list(mean.intercept=...)} and \code{control.fixed=list(mean=...)}. Same as \code{mean.intercept.loc}.
#' @param prec.loc the prior precision for all the Gaussian additive term for each node. INLA argument \code{control.fixed=list(prec.intercept=...)} and \code{control.fixed=list(prec=...)}. Same as \code{prec.intercept.loc}.
#' @param loggam.shape.loc the shape parameter in the Gamma distribution prior for the precision in a Gaussian node. INLA argument \code{control.family=list(hyper = list(prec = list(prior="loggamma",param=c(loggam.shape, loggam.inv.scale))))}.
#' @param loggam.inv.scale.loc the inverse scale parameter in the Gamma distribution prior for the precision in a Gaussian node. INLA argument \code{control.family=list(hyper = list(prec = list(prior="loggamma",param=c(loggam.shape, loggam.inv.scale))))}.
#' @param verbose.loc FALSE to not print additional output.
#' @param nthreads number of threads to use for INLA. Default is \code{fit.control[["ncores"]]} or \code{build.control[["ncores"]]} which is the number of cores specified in \code{control} and defaults to 1.
#'
#' @return If INLA failed, FALSE or an error is returned. Otherwise, the direct output from INLA is returned.
#' @family Bayes
#' @keywords internal
calc.node.inla.glmm <-
  function(child.loc = NULL,
           dag.m.loc = NULL,
           data.df.loc = NULL,
           data.dists.loc = NULL,
           ntrials.loc = NULL,
           exposure.loc = NULL,
           compute.fixed.loc = NULL,
           mean.intercept.loc = NULL,
           prec.intercept.loc = NULL,
           mean.loc = NULL,
           prec.loc = NULL,
           loggam.shape.loc = NULL,
           loggam.inv.scale.loc = NULL,
           verbose.loc = FALSE,
           nthreads = NULL) {

    if (nthreads == 1) {
      INLA::inla.setOption("num.threads", "1:1")
    } else if (nthreads > 1) {
      # inlathreads <- paste0(nthreads, ":1")
      # INLA::inla.setOption("num.threads", inlathreads)
      if (verbose.loc) {
        message("Nested parallelism detected. Limiting INLA (inner loop) to 1 thread to prevent unexpected behaviour.\n")
      }
      INLA::inla.setOption("num.threads", "1:1")
    } else if (nthreads < 1) {
      stop("invalid number of threads for INLA")
    }
    if (verbose.loc) {message("INLA threads (outer:inner) set to ", INLA::inla.getOption("num.threads"), "\n")}

    #print(data.df.loc);
    #print(group.var);
    group.var <-
      names(data.df.loc)[length(names(data.df.loc))]
    ## group variable is always the last column
    ## 1. get the formula part of the call - create a string of this
    if (length(which(dag.m.loc[child.loc,-child.loc] == 1)) == 0) {
      ## independent node
      str.eqn.str <-
        paste(colnames(dag.m.loc)[child.loc], "~1+")

    } else {
      ## have some covariate

      if (dim(dag.m.loc)[1] == 2) {
        ## special case - 2x2 DAG and so names are not retained when -child.loc
        str.eqn.str <-
          paste(colnames(dag.m.loc)[child.loc],
                "~",
                colnames(dag.m.loc)[-child.loc],
                "+",
                sep = "")

      } else {
        str.eqn.str <-
          paste(colnames(dag.m.loc)[child.loc],
                "~",
                paste(names(which(
                  dag.m.loc[child.loc,-child.loc] == 1
                )), collapse = "+", sep = ""),
                "+",
                sep = "")
      }
    }
    #cat(str.eqn.str,"\n");
    #stop("");
    ## new part - add in the function for latent variables
    str.eqn.str <-
      paste(
        str.eqn.str,
        "f(",
        group.var,
        ",model=\"iid\",hyper=list(theta=list(prior=\"loggamma\",param=c(",
        loggam.shape.loc,
        ",",
        loggam.inv.scale.loc,
        ")))),\n",
        sep = ""
      )

    #print(str.eqn.str);
    #stop("");
    ## 2. data set
    str.data <- "data=data.df.loc,"


    ## 3. family
    str.family <-
      paste("family=\"", data.dists.loc[[child.loc]], "\",", sep = "")


    ## 4. additional parameter for number of trials (binomial) or exposure (poisson)
    str.extra <- ""

    if (data.dists.loc[[child.loc]] == "binomial") {
      str.extra <- paste("Ntrials=ntrials.loc,", sep = "")
    }
    if (data.dists.loc[[child.loc]] == "poisson") {
      str.extra <- paste("E=exposure.loc,", sep = "")
    }
    if (data.dists.loc[[child.loc]] == "gaussian") {
      str.extra <-
        paste(
          "control.family=list(hyper = list(prec = list(prior=\"loggamma\",param=c(",
          loggam.shape.loc,
          ",",
          loggam.inv.scale.loc,
          ")))),\n",
          sep = ""
        )
    }

    ## 5. get the full command
    res <- NULL

    start.str <- "res <- INLA::inla("

    end.str <-
      paste(
        "\ncontrol.fixed=list(mean.intercept=",
        mean.intercept.loc,
        ",\n",
        "prec.intercept=",
        prec.intercept.loc,
        ",\n",
        "mean=",
        mean.loc,
        ",\n",
        "prec=",
        prec.loc,
        ",\n",
        "compute=",
        compute.fixed.loc,
        "))\n",
        sep = ""
      )


    r <- NULL

    full.command <-
      paste(
        "r <- try(",
        start.str,
        str.eqn.str,
        str.data,
        str.family,
        str.extra,
        end.str,
        ",silent=TRUE)",
        sep = ""
      )


    ## 6. some debugging - if requested
    if (verbose.loc) {
      cat("commands which are parsed and sent to inla().\n")

      print(full.command)
    }

    ## 7. now run the actual command - parse and eval - is parsed in current scope and so data.df exists here
    eval(parse(text = full.command))

    if (inherits(r, what = "try-error")) {
      warning(r)
      return(FALSE)
    } else if (length(r) == 1) {
      ### INLA failed
      warning("INLA failed\n")
      return(FALSE)
    } else {
      ## 8. return the results
      if (!compute.fixed.loc) {
        ## only want marginal likelihood
        return(res$mlik[1])
        ## n.b. [1] is so we choose the integrated rather than Gaussian version - debateable which to choose

      } else {
        ## alternatively get *all* the output from inla() - copious
        return(res)
      }
    }

  } ## end of function
