# Script that runs code extracted from the buildScoreCache() function
# It runs all checks and data cleaning before calling the buildScoreCache.mle()

setupScoreCache <- function(data.df = NULL,
                            data.dists = NULL,
                            method = "bayes",
                            group.var = NULL,
                            adj.vars = NULL,
                            cor.vars = NULL,
                            dag.banned = NULL,
                            dag.retained = NULL,
                            max.parents = NULL,
                            which.nodes = NULL,
                            defn.res = NULL,
                            centre = TRUE,
                            dry.run = FALSE,
                            control = NULL,
                            verbose = FALSE,
                            debugging = FALSE){

  ## start tests
  # Check verbose
  if(!any(verbose %in% c(TRUE, FALSE))){
    stop(paste("'verbose' is not provided but should be TRUE/FALSE."))
  }

  # Check data
  if(!is.null(data.df)){
    mylist <- check.valid.data(data.df = data.df, data.dists = data.dists, group.var = group.var)
  } else if (is.na(data.df)){
    stop("'data.df' is NA but must be provided.")
  } else {
    stop(paste("'data.df' is not provided."))
  }

  # Check dists
  if(!is.null(data.dists)){
    data.dists <- validate_dists(data.dists = data.dists, returnDists=TRUE)
  } else if (is.na(data.dists)){
    stop("'data.dists' is NA but must be provided.")
  } else {
    stop(paste("'data.dists' is not provided."))
  }

  # Check method
  if(!is.null(method)){
    method <- tolower(method)
    if (is.na(method)){
      stop("'method' is NA but must be provided.")
    } else if(!(method %in% c("bayes", "mle"))){
      stop("`method` is unknown.")
    }
  } else {
    stop(paste("'method' is not provided."))
  }

  # Check group.var
  val_groups <- check.valid.groups(group.var = group.var, data.df = data.df, cor.vars = cor.vars, verbose = verbose)

  # Check cor.vars
  val_corvars <- check.valid.groups(group.var = group.var, data.df = data.df, cor.vars = cor.vars, verbose = verbose)

  if (!is.null(group.var) && !is.null(cor.vars) && !is.null(which.nodes)){
    if (group.var %in% colnames(data.df)[which.nodes]){
      stop("`group.var` should not be among `which.nodes`.")
    } else {
      data.df <- val_groups[["data.df"]]
    }
  } else {
    #TODO: think about this situation...
  }

  # Check adj.vars
  if((!is.null(adj.vars) & !is.null(cor.vars)) & !(is.null(cor.vars[adj.vars]))){stop("cor.vars contains adj.vars, please remove them")}

  # Check dag.banned
  dag.banned <- check.valid.dag(dag = dag.banned, data.df = data.df, is.ban.matrix = TRUE, group.var = group.var)

  # Check dag.retained
  dag.retained <- check.valid.dag(dag = dag.retained, data.df = data.df, is.ban.matrix = FALSE, group.var = group.var) # is.ban.matrix = FALSE to check for cycles

  ###
  # Check max.parents
  # which.nodes: Limits the number of variables. Selects the variables for consideration. Integer vector of column indeces.
  # max.parents: Limits the number of allowed parent variables. single, Integer, list or NULL. Output: integer vector of max.parents corresponding to the column index.
  # defn.res: User specified child-parent combinations. Should not disagree with which.nodes and the max.number of parents.

  # make here max.parents a vector first? This would disentangle the handling of the different input types of max.parents?
  max.parents.orig <- max.parents # store original max.parents for later
  if(is.list(max.parents)){
    # make numeric vector out of list.
    max.parents <- unlist(max.parents, use.names = FALSE)

    # check if all elements were unlisted
    if(length(max.parents) != length(max.parents.orig)){
      stop("Something is wrong in the max.parents list. Any unallowed NA or NULL?")
    }
  }

  # max.parent is NULL or numeric (vector) from here on.
  if (!is.null(defn.res)) {
    if (!is.null(which.nodes)) {
      if (!is.null(max.parents)) {
        # check max.parents on trimmed data set and trimmed distribution
        max.parents <- check.valid.parents(data.df = data.df[which.nodes], max.parents = max.parents, group.var = group.var)
      } else {
        # no max.parents provided
        max.p <- length(data.dists[which.nodes]) # not (n-1) here, to raise warning downstream.

        # check max.p
        max.parents <- check.valid.parents(data.df = data.df[which.nodes], max.parents = max.p, group.var = group.var)
        if(verbose){message("`max.parents` was provided as NULL. Therefore, I assume max.parents equal to the number of possible nodes.")}
      }
      # Check consistency of which.nodes and max.parents
      if(length(which.nodes) < max(unique(max.parents))){
        stop("`which.nodes` selected less nodes than max.parents assumes.")
      }
    } else {
      # no which.nodes provided
      if (!is.null(max.parents)) {
        # check max.parents
        max.parents <- check.valid.parents(data.df = data.df, max.parents = max.parents, group.var = group.var) # this works if max.parents is a single numeric, a list or a numeric vector.
      } else {
        # no max.parents provided
        max.p <- length(data.dists) # not (n-1) here, to raise warning downstream.

        # check max.p
        max.parents <- check.valid.parents(data.df = data.df, max.parents = max.p, group.var = group.var)
        if(verbose){message("`max.parents` was provided as NULL. Therefore, I assume max.parents equal to the number of possible nodes.")}
      }
    }
    # check if defn.res agrees with which.nodes
    which.nodes.defnres <- unique(defn.res$children)
    if (which.nodes.defnres != which.nodes){
      stop("`defn.res` and `which.nodes` provided but can only use either or.")
    }
    # check if defn.res agrees with max.parents
    max.parents.defnres <- max(apply(defn.res[["node.defn"]],1,sum))
    if (max.parents.defnres != max.parents){
      stop("`defn.res` and `max.parents` provided but can only use one of each.")
    }
  } else {
    # no defn.res provided
    if (!is.null(which.nodes)) {
      if (!is.null(max.parents)) {
        # check max.parents on trimmed data set and trimmed distribution
        max.parents <- check.valid.parents(data.df = data.df[which.nodes], max.parents = max.parents, group.var = group.var)
      } else {
        # no max.parents provided
        max.p <- length(data.dists[which.nodes]) # not (n-1) here, to raise warning downstream.

        # check max.p
        max.parents <- check.valid.parents(data.df = data.df[which.nodes], max.parents = max.p, group.var = group.var)
        if(verbose){message("`max.parents` was provided as NULL. Therefore, I assume max.parents equal to the number of possible nodes.")}
      }
      # Check consistency of which.nodes and max.parents
      if(length(which.nodes) < max(unique(max.parents))){
        stop("`which.nodes` selected less nodes than max.parents assumes.")
      }
    } else {
      # no which.nodes provided
      if (!is.null(max.parents)) {
        # check max.parents
        max.parents <- check.valid.parents(data.df = data.df, max.parents = max.parents, group.var = group.var) # this works if max.parents is a single numeric, a list or a numeric vector.
      } else {
        # no max.parents provided
        max.p <- length(data.dists) # not (n-1) here, to raise warning downstream.

        # check max.p
        max.parents <- check.valid.parents(data.df = data.df, max.parents = max.p, group.var = group.var)
        if(verbose){message("`max.parents` was provided as NULL. Therefore, I assume max.parents equal to the number of possible nodes.")}
      }
    }
  } # EOF checking max.parents, which.nodes, defn.res
  # This outputs max.parents as numeric (vector) with a maximal value for max.parents equal to the number of nodes in data.df or data.df[which.nodes] respectively (max.parents <= n).

  # If max.parents = max number of nodes, correct for (n-1)
  if(!is.null(which.nodes)){
    nvars <- length(data.dists[which.nodes])
  } else if(is.null(which.nodes)){
    nvars <- length(data.dists)
  } else {
    stop("`which.nodes` is invalid.")
  }

  if (length(max.parents) == 1){
    # max.parents is single integer
    if(max.parents>=nvars){
      max.parents <- nvars-1
      warning(paste("`max.parents` >= no. of variables. I set it to (no. of variables - 1)=", max.parents))
    }
  } else {
    # max.parents is numeric vector
    counter <- 0
    for (i in 1:length(max.parents)){
      if (max.parents[i] >= nvars) {
        max.parents[i] <- max.parents[i]-1
        counter <- counter+1
      }
    }
    if(counter >0){warning("Some values of `max.parents` >= no. of variables. I set them to (no. of variables - 1): ", max.parents)}
  }

  # Check which.nodes and defn.res
  if (is.null(defn.res)){
    if (is.null(which.nodes)){
      # make which.nodes = all.nodes
      which.nodes <- check.which.valid.nodes(data.df = data.df, which.nodes = which.nodes, group.var = group.var)
    } else {
      # check the provided which.nodes
      which.nodes <- check.which.valid.nodes(data.df = data.df, which.nodes = which.nodes, group.var = group.var)
    }
  } else {
    if (is.null(which.nodes)){
      # If user supplied children and parent combinations, make which.nodes based on them.
      which.nodes <- unique(defn.res$children)
    } else {
      stop("`defn.res` and `which.nodes` provided but can only use either or.")
    }
  }

  # Check centre
  if(!any(centre %in% c(TRUE, FALSE))){
    if(is.null(centre)){
      centre <- TRUE
      warning("`centre` is not provided. I set it to TRUE.")
    } else {
      stop(paste("`centre` is not provided but should be TRUE/FALSE."))
    }
  } else {
    # centre is provided as TRUE/FALSE
    # Check if there are nodes to centre or if there are no gaussian nodes
    if(!any(data.dists == "gaussian")) {
      centre <- FALSE
      if(verbose){message("No gaussian nodes to centre. I set `centre` to FALSE.")}
    } else {
      if(verbose){message("I set `centre` to: ", centre)}
      centre <- centre
    }
  }

  # Check dry.run
  if(!any(dry.run %in% c(TRUE, FALSE))){
    stop(paste("'dry.run' is not provided but should be TRUE/FALSE."))
  }

  # Check control args
  # if any arg matches a possible control arg from build.control(), warn and use it.
  build.control.args <- names(formals(build.control))[-which(names(formals(build.control))=="method")] # allow method to be provided as it has a different meaning here.
  provided.args <- names(match.call()[-1]) # remove the first element of match.call() which is empty.
  if(any(provided.args %in% build.control.args)){
    warning(paste("Some arguments match possible control arguments from `build.control()`.
                  I will use the provided arguments. Please use `control=build.control(...)` instead in the future."))
    ambiguous.args <- provided.args[which(provided.args %in% build.control.args)]
    for (i in 1:length(ambiguous.args)){
      control[[ambiguous.args[i]]] <- match.call()[-1][[ambiguous.args[i]]]
    }
  }
  ctrl <- check.valid.buildControls(control = control, method = method, verbose = verbose)

  if("max.mode.error" %in% names(ctrl)){
    if(ctrl[["max.mode.error"]]==0) {
      force.method <- "C"
    } else if(ctrl[["max.mode.error"]]==100) {
      force.method <- "INLA"
    } else {
      force.method <- "notset"
    }
  } else {
    # no max.mode.error among control list
    force.method <- "notset"
  }

  ## Check consistency of dag.retain and dag.banned
  # check retain does not ask for more arcs to be retained than allowed in max.parents
  # number of parents per node to retain
  max.retain <- apply(dag.retained,1,sum)
  if(length(which( (max.retain>max(unique(max.parents))) == TRUE))>0){stop("'dag.retained' is inconsistent with max.parents!")}

  # check that arcs than are banned are also not retained
  if(length(which(which(as.integer(dag.banned)==1)%in%which(as.integer(dag.retained)==1)==TRUE))>0){stop("'dag.banned' and 'dag.retained' are inconsistent!")}

  # returns ammended data.df and suitable variables
  list.group.var <- check.valid.groups(group.var = group.var, data.df = data.df, cor.vars = cor.vars)
  ## int vect of variables to be treated as grouped indexed from 1
  grouped.vars <- list.group.var$grouped.vars
  ## int vector of group membership ids
  group.ids <- list.group.var$group.ids
  ## this has removed the grouping variable from data.df
  data.df <- list.group.var$data.df

  
  return(list(data.df = data.df, data.dists = data.dists,
	      max.parents = max.parents, adj.vars = adj.vars,
	      cor.vars = cor.vars, dag.banned = dag.banned,
	      dag.retained = dag.retained, which.nodes = which.nodes,
	      centre = centre, defn.res = defn.res, dry.run = dry.run,
	      verbose = verbose, debugging = debugging,
	      force.method = force.method, group.var = group.var,
	      grouped.vars = grouped.vars, group.ids = group.ids, control = ctrl))
}
