#' Perform inference with ABN
#'
#' Main function to predict the distribution of a node in a fitted ABN graph
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param dag An adjacency matrix (can be the output of the function mostProbable()).
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param hypothesis Node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param plot TRUE/FALSE to indicate if the predicted distribution has to be plotted.
#' @return A list containing the predicted distribution of the hypothesis and predicted distributions of the upstream nodes.
#' @import igraph
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#'
#' hypothesis <- "g2"
#'
#' evidence <- list("p1" = 3)
#'
#' predictions <- predictABN(data, mydists, dag, myfit, hypothesis, evidence)
#' str(predictions)
#' predictions$prediction_hypothesis # the predicted distribution of g2
#' @export
#'
predictABN <- function(data, dists, dag, fit, hypothesis = NULL, evidence = NULL, plot = FALSE){

  # some checks
  if (ncol(data) != length(dists) || ncol(data) != ncol(dag) || length(dists) != ncol(dag)){
    stop("The number of nodes/variables do not correspond.")
  }
  if (!all(colnames(data) %in% names(dists)) || !all(colnames(data) %in% colnames(dag)) || !all(names(dists) %in% colnames(dag))){
    stop("The names of the nodes/variables in data, dists and dag do not correspond.")
  }
  graph <- graph_from_adjacency_matrix(t(dag))
  node_order <- names(topo_sort(graph, mode="out"))

  # Step 0: check the evidence and the data
  check_data(data, dists,fit)
  evidence <- check_evidence(data, dists, hypothesis, evidence)
  if (!is.null(hypothesis)){
    mb <- find_MB(graph, hypothesis)
  } else {
    mb <- NULL
  }
  if (all(mb %in% names(evidence)) && !is.null(hypothesis)){
    message("The Markov blanket of the node to infer is enterily known.")
    prediction <- predict_node_from_parent(data, dists, graph, fit, node = hypothesis, evidence)
    predictions <- list(prediction)
    names(predictions) <- hypothesis
    prediction_updated <- predict_node_from_children(data, dists, graph, fit, node = hypothesis, evidence, predictions = predictions)
    return(list(prediction_hypothesis = prediction_updated, predictions = NULL))
  } else {
    # Step 1: from top to bottom
    predictions <- list()
    predictions_names <- c()
    for (i in (1:length(node_order))){
      prediction <- predict_node_from_parent(data, dists, graph, fit, node = node_order[i], evidence, predictions)
      predictions <- c(predictions, list(prediction))
      predictions_names <- c(predictions_names, node_order[i])
      names(predictions) <- predictions_names
    }

    # Step 2: from bottom to top
    if (!is.null(hypothesis)){
      node_max <- which(node_order == hypothesis)
    } else {
      node_max <- 1
    }
    for (i in (length(node_order):node_max)){
      prediction <- predict_node_from_children(data, dists, graph, fit, node = node_order[i], evidence, predictions)
      predictions[[i]] <- prediction
    }

    # Step 3: Plot the posterior distribution
    if (plot==TRUE){
      g <- plotPosteriorDistrib(predictions,hypothesis,dists)
      print(g)
    }
    if (!is.null(hypothesis)){
      results <- list(prediction_hypothesis = predictions[[hypothesis]], predictions = predictions)
    } else {
      results <- list(prediction_hypothesis = NULL, predictions = predictions)
    }
    return(results)
  }
}

#' Check the evidence
#'
#' This function checks the format of the evidence
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param hypothesis Node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @return A list containing the evidence in the right format
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' hypothesis <- "g2"
#' evidence <- list("b2" = "y", b5 = "n")
#'
#' evidence <- check_evidence(data, mydists, hypothesis, evidence)
#'
#' @export
check_evidence <- function(data, dists, hypothesis, evidence){
  if (length(evidence)>0){
    # at least one evidence
    if (hypothesis %in% names(evidence)){
      print(paste0("Node ",hypothesis," is already known (evidence)."))
      evidence <- evidence[hypothesis]
    } else {
      evidence.to.remove <- c()
      for (i in (1:length(evidence))){
        dist.evidence <- dists[[names(evidence)[i]]]
        if (dist.evidence == "binomial" || dist.evidence == "multinomial"){
          if (! evidence[[i]] %in% levels(data[[names(evidence[i])]])){
            warning(paste0("Evidence ",names(evidence[i])," does not have an expected value. It should be either ",paste(levels(data[[names(evidence[i])]]),collapse=", "),". It will be discarded."))
            evidence.to.remove <- c(evidence.to.remove,names(evidence[i]))
          }
        } else if (dist.evidence == "poisson" || dist.evidence == "gaussian"){
          if (is.character(evidence[[i]])){
            warning(paste0("Evidence ",names(evidence[i])," is a string but should be numeric. It will be discarded."))
            evidence.to.remove <- c(evidence.to.remove,names(evidence[i]))
          }
        }
      }
      if (length(evidence.to.remove)>0){
        evidence <- evidence[setdiff(names(evidence),evidence.to.remove)]
      }

      # rewrite the binomial evidences
      if (length(which(dists[names(evidence)] %in% c("binomial","multinomial")))>0){
        # at least one evidence is a bi/multinomial node
        for (i in (1:length(which(dists[names(evidence)]%in% c("binomial","multinomial"))))){
          node <- names(evidence)[which(dists[names(evidence)] %in% c("binomial","multinomial"))][i]
          if (is.character(evidence[[node]])){
            # transform it to factor
            evidence[[node]] <- factor(evidence[[node]],levels=levels(data[[node]]))
          }
        }
      }
    }
  }
  return(evidence)
}

check_data <- function(data, dists,fit){
  data.bin <- data %>% dplyr::select(names(dists)[grep("b",names(dists))])

  level.length <- sapply(data.bin, function(b){
    length(levels(b))
  })

  if (length(which(level.length==1))>0){
    stop(paste0("Binomial node ",names(level.length)[which(level.length==1)]," does not have the right number of levels (2). Consider adding one level (data$bin.node <- factor(data$bin.node,levels=c(0,1)) before running the code."))
  }

  data.multi <- data %>% dplyr::select(names(dists)[grep("m",names(dists))])

  level.length.multi <- sapply(data.multi, function(m){
    length(levels(m))
  })

#  true.levels.multi <- sapply(names(level.length.multi),function(m){
#    length(grep("intercept",names(fit[[m]])))+1
#  })

#  multi.NA <- names(level.length.multi)[which(is.na(match(unlist(level.length.multi),unlist(true.levels.multi))))]
#  if (length(multi.NA)>0){
  if ((length(which(level.length.multi<3))>0)){
    #stop(paste0("Multinomial node ",multi.NA," does not have the expected number of levels (",true.levels.multi[multi.NA],"). Consider adding one level (data$multi.node <- factor(data$multi.node,levels=c(0,1)) before running the code."))
    stop(paste0("Multinomial node ",names(level.length)[which(level.length.multi<3)]," does not have the expected number of levels. Consider adding one or more levels (data$multi.node <- factor(data$multi.node,levels=c(0,1,2)) before running the code."))
  }
}

#' Find the Markov Blanket of a node
#'
#' Identifies the Markov Blanket of a specific node within the DAG.
#' The Markov Blanket includes the node's parents, its children, and any other
#' parents of those children.
#'
#' @param graph An igraph object representing the DAG.
#' @param hypothesis A character string specifying the name of the node for which to find the Markov Blanket.
#'
#' @return A character vector containing the names of the nodes in the Markov Blanket.
#'
#' @export
find_MB <- function(graph, hypothesis){
  parents <- find_parents(node = hypothesis,graph = graph)
  children <- find_children(node = hypothesis, graph=graph)
  parents_children <- sapply(children, function(child){
    find_parents(node = child,graph = graph)
  })
  MB <- setdiff(unique(c(parents, children, unlist(parents_children))),hypothesis)
  return(MB)
}

#' Perform upstream inference with ABN
#'
#' Main function to predict the distribution of a node given its parents only
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param graph A dag (igraph object).
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param predictions The estimated predictions of the upstream nodes (can be empty if the parents of the node to predict are evidence).
#' @return The predicted distribution of the node of interest.
#' @import igraph
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#' graph <- igraph::graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#'
#' node <- "g2"
#'
#' evidence <- list("p1" = 3, "g1" = 5, "b2" = "y")
#'
#' predictions <- predict_node_from_parent(data, mydists, graph, myfit, node, evidence)
#' predictions  # the predicted distribution of g2
#' @export
#'
predict_node_from_parent <- function(data, dists, graph, fit, node, evidence, predictions = NULL){
  parents <- find_parents(graph, node)
  if (is.null(predictions)){
    if (!all(parents %in% names(evidence))){
      # all parents are not evidence
      stop("Not enough information about the upstream nodes.")
    }
    predictions <- evidence
  }
  results <- switch(dists[[node]],
                    "poisson"  = predict_node_from_parent_poisson(data, dists, fit, node, evidence, parents, predictions),
                    "gaussian" = predict_node_from_parent_gaussian(data, dists, fit, node, evidence, parents, predictions),
                    "binomial" = predict_node_from_parent_binomial(data, dists, fit, node, evidence, parents, predictions),
                    "multinomial" = predict_node_from_parent_multinomial(data, dists, fit, node, evidence, parents, predictions),
                    stop("Unknown node distribution type: ", dists[[node]]) # Error handling
  )
  return(results)
}

#' Perform downstream inference with ABN
#'
#' Main function to predict the distribution of a node given its children and the children's parents
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param graph A dag (igraph object).
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param predictions The estimated predictions of the downstream nodes (must contain at least a first prediction of the node to predict if the children and the children's parents are evidence).
#' @import igraph
#' @return The predicted distribution of the node of interest.
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#' graph <- igraph::graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#'
#' node <- "g2"
#'
#' evidence <- list("g1" = 3, "b5" = "y")
#'
#' predictions <- list("g2" = c(0,1)) # a first estimate of the node g2
#'
#' predictions <- predict_node_from_children(data, mydists, graph, myfit, node, evidence, predictions)
#' predictions  # the predicted distribution of g2
#' @export
#'
predict_node_from_children <- function(data, dists, graph, fit, node, evidence, predictions){
  if (!node %in% names(predictions)){
    stop("Predictions must contain at least a first prediction of the node to predict.")
  }

  children <- find_children(graph, node)
  if (length(children)==0){
    results <- predictions[[node]]
    return(results)
  }

  parents_children <- unique(unlist(sapply(children,function(l){
    find_parents(graph, node = l)
  })))

  if (!all(c(children,parents_children) %in% names(predictions))){
    nodes <- setdiff(c(children,parents_children),names(predictions))
    if (!all(nodes %in% names(evidence))){
      stop("Not enough information about the downstream nodes.")
    }
    predictions <- c(predictions,evidence)
  }

  Results <- lapply(children, function(child) {
    parents <- find_parents(graph, child)
    dist_type <- dists[[child]]

    switch(dist_type,
           "gaussian"  = predict_node_from_children_gaussian(data, dists, fit, node, evidence, child, parents, predictions),
           "poisson"   = predict_node_from_children_poisson(data, dists, fit, node, evidence, child, parents, predictions),
           "binomial"  = predict_node_from_children_binomial(data, dists, fit, node, evidence, child, parents, predictions),
           "multinomial" = predict_node_from_children_multinomial(data, dists, fit, node, evidence, child, parents, predictions),
           stop(paste("Unsupported distribution type:", dist_type))
    )
  })

  raw_matrix <- do.call(rbind, Results)

  switch(dists[[node]],
         "gaussian" = {
           res_mean <- mean(raw_matrix[, 1])
           res_var  <- mean(raw_matrix[, 2]) / length(children)
           return(c(res_mean, res_var))
         },
         "binomial" = {
           node_levels <- levels(data[[node]])

           if (is.null(ncol(raw_matrix)) || ncol(raw_matrix) == 1) {
             res <- numeric(length(node_levels))
             names(res) <- node_levels

             chosen_index <- as.integer(raw_matrix[1, 1])
             res[chosen_index] <- 1
             return(res)
           }

           p_level2 <- mean(raw_matrix[, 2])
           res <- c(1 - p_level2, p_level2)
           names(res) <- node_levels
           return(res)
         },
         "multinomial" = {
           node_levels <- levels(data[[node]])

           if (is.null(ncol(raw_matrix)) || ncol(raw_matrix) == 1) {
             res <- numeric(length(node_levels))
             names(res) <- node_levels

             chosen_index <- as.integer(raw_matrix[1, 1])
             res[chosen_index] <- 1
             return(res)
           }
           res <- colMeans(raw_matrix)
           res <- res / sum(res)
           names(res) <- node_levels
           return(res)
         },
         "poisson" = {
           return(mean(raw_matrix[, 1]))
         },
         # Default case for aggregation
         stop(paste("Aggregation logic missing for distribution:", dists[[node]]))
  )
}

#' Perform upstream inference with ABN
#'
#' Main function to predict the distribution of a Poisson node given its parents only
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param parents The parents of the node to predict.
#' @param predictions The estimated predictions of the upstream nodes  (can be empty if the parents of the node to predict are evidence).
#' @return The predicted distribution of the node of interest.
#' @import igraph
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#' graph <- igraph::graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#'
#' node <- "p2"
#' parents <- find_parents(graph,node)
#'
#' evidence <- list("b1" = "y", "p1" = 3)
#'
#' predictions <- predict_node_from_parent_poisson(data, mydists, myfit, node, evidence, parents)
#' predictions  # the predicted distribution of p1
#' @export
#'
predict_node_from_parent_poisson <- function(data, dists, fit, node, evidence, parents, predictions = NULL){
  if (dists[[node]] != "poisson"){
    stop("The node to predict should follow a Poisson distribution.")
  }

  if (node %in% names(evidence)){
    # node is an evidence
    node_hat <- evidence[[node]]
  } else {
    if (length(parents)==0){
      node_hat <- predict_root(data, dists, node)
    } else {
      eq <- fit[[node]]
      eq_names <- colnames(eq) %||% names(eq)
      eq <- setNames(as.vector(eq), gsub(".*\\|", "", eq_names))

      bin.nodes <- intersect(names(dists)[which(dists %in% c("binomial","multinomial"))],parents)
      if (length(bin.nodes)>0){
        other.nodes <- parents[-which(parents %in% bin.nodes)]
      } else {
        other.nodes <- parents
      }

      predictions_tmp <- predictions[other.nodes]
      predictions_tmp <- lapply(predictions_tmp, function(l){ l[1] })
      continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
      names(continuous_part) <- c()

      if (length(bin.nodes)>0){
        probabilities <- predictions[bin.nodes]

        bin.nodes.evidence <- intersect(names(evidence),bin.nodes)
        if (length(bin.nodes.evidence)>0){
          for (i in (1:length(bin.nodes.evidence))){
            proba_tmp <- rep(0,length(levels(probabilities[[bin.nodes.evidence[i]]])))
            names(proba_tmp) <- levels(data[[bin.nodes.evidence[i]]])
            proba_tmp[grep(probabilities[[bin.nodes.evidence[i]]],names(proba_tmp))] <- 1
            probabilities[[bin.nodes.evidence[i]]] <- proba_tmp
          }
        }
        names(probabilities) <- bin.nodes
        levels_list <- lapply(probabilities,function(x) names(x))
        combinations <- expand.grid(levels_list)

        bin_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "binomial")]
        multi_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "multinomial")]

        if(length(bin_parents) > 0) {
          mat_bin <- as.matrix(sapply(combinations[bin_parents], function(x) as.numeric(as.character(x))))
        } else {
          mat_bin <- NULL
        }

        if(length(multi_parents) > 0) {
          df_multi <- combinations[multi_parents]
          for(n in multi_parents) df_multi[[n]] <- factor(df_multi[[n]], levels = names(probabilities[[n]]))
          mat_multi <- model.matrix(~ ., data = df_multi)
          mat_multi <- mat_multi[, colnames(mat_multi) != "(Intercept)", drop = FALSE]
        } else {
          mat_multi <- NULL
        }

        dummy_matrix <- cbind(mat_bin, mat_multi)

        combinations_tmp <- dummy_matrix %*% eq[colnames(dummy_matrix)]

        proba_cond_values <- exp(continuous_part + rowSums(combinations_tmp))

        prob_grid <- expand.grid(probabilities)
        combination_probabilities <- apply(prob_grid, 1, prod)

        node_hat <- sum(proba_cond_values * combination_probabilities)
      } else {
        node_hat <- exp(continuous_part)
      }
    }
  }
  return(node_hat = node_hat)
}

#' Perform upstream inference with ABN
#'
#' Main function to predict the distribution of a Gaussian node given its parents only
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param parents The parents of the node to predict.
#' @param predictions The estimated predictions of the upstream nodes  (can be empty if the parents of the node to predict are evidence).
#' @return The predicted distribution of the node of interest.
#' @import igraph
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#' graph <- igraph::graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#'
#' node <- "g2"
#' parents <- find_parents(graph,node)
#'
#' evidence <- list("p1" = 1, "g1" = 2, "b2" = "y")
#'
#' predictions <- predict_node_from_parent_gaussian(data, mydists, myfit, node, evidence, parents)
#' predictions  # the predicted distribution of g2
#' @export
#'
predict_node_from_parent_gaussian <- function(data, dists, fit, node, evidence, parents, predictions = NULL){
  if (dists[[node]] != "gaussian"){
    stop("The node to predict should follow a Gaussian distribution.")
  }
  if (node %in% names(evidence)){
    # node is an evidence
    node_hat <- c(evidence[[node]],var(data[[node]]))
  } else {
    if (length(parents)==0){
      # no parents
      node_hat <- predict_root(data, dists, node)
      node_hat <- c(node_hat,fit[[node]][2])
    } else {
      eq <- fit[[node]]
      eq_names <- colnames(eq) %||% names(eq)
      eq <- setNames(as.vector(eq), gsub(".*\\|", "", eq_names))

      node_sigma_sq <- 1 / eq["precision"]
      bin.nodes <- intersect(names(dists)[which(dists %in% c("binomial","multinomial"))],parents)
     if (length(bin.nodes)>0){
        other.nodes <- parents[-which(parents %in% bin.nodes)]
      } else {
        other.nodes <- parents
      }

      predictions_tmp <- predictions[other.nodes]
      predictions_tmp <- lapply(predictions_tmp,function(l){
        l[1]
      })
      continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
      names(continuous_part) <- c()

      if (length(bin.nodes)>0){
        probabilities <- predictions[bin.nodes]

        bin.nodes.evidence <- intersect(names(evidence),bin.nodes)

        if (length(bin.nodes.evidence)>0){
          # at least one bin nodes is an evidence
          for (i in (1:length(bin.nodes.evidence))){
            proba_tmp <- rep(0,length(levels(probabilities[[bin.nodes.evidence[i]]])))
            names(proba_tmp) <- levels(data[[bin.nodes.evidence[i]]])
            proba_tmp[grep(probabilities[[bin.nodes.evidence[i]]],names(proba_tmp))] <- 1
            probabilities[[bin.nodes.evidence[i]]] <- proba_tmp
          }
        }
        names(probabilities) <- bin.nodes
        levels_list <- lapply(probabilities,function(x) names(x))
        combinations <- expand.grid(levels_list)

        bin_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "binomial")]
        multi_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "multinomial")]

        if(length(bin_parents) > 0) {
          mat_bin <- as.matrix(sapply(combinations[bin_parents], function(x) as.numeric(as.character(x))))
        } else {
          mat_bin <- NULL
        }

        if(length(multi_parents) > 0) {
          df_multi <- combinations[multi_parents]
          for(n in multi_parents) df_multi[[n]] <- factor(df_multi[[n]], levels = names(probabilities[[n]]))
          mat_multi <- model.matrix(~ ., data = df_multi)
          mat_multi <- mat_multi[, colnames(mat_multi) != "(Intercept)", drop = FALSE]
        } else {
          mat_multi <- NULL
        }

        dummy_matrix <- cbind(mat_bin, mat_multi)

        combinations_tmp <- dummy_matrix %*% eq[colnames(dummy_matrix)]

        proba_cond_values <- continuous_part + as.vector(combinations_tmp)

        prob_grid <- expand.grid(probabilities)
        combination_probabilities <- apply(prob_grid, 1, prod)

        node_hat <- c(sum(proba_cond_values * combination_probabilities),node_sigma_sq)
      } else {
        node_hat <- c(continuous_part,node_sigma_sq)
      }
    }
  }
  return(node_hat = node_hat)
}

#' Perform upstream inference with ABN
#'
#' Main function to predict the distribution of a binomial node given its parents only
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param parents The parents of the node to predict.
#' @param predictions The estimated predictions of the upstream nodes  (can be empty if the parents of the node to predict are evidence).
#' @return The predicted distribution of the node of interest.
#' @import igraph
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#' graph <- igraph::graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#'
#' node <- "b3"
#' parents <- find_parents(graph,node)
#'
#' evidence <- list("g1" = 4, "b1" = "y", "b2" = "y")
#'
#' predictions <- predict_node_from_parent_binomial(data, mydists, myfit, node, evidence, parents)
#' predictions  # the predicted distribution of b2
#' @export
#'
predict_node_from_parent_binomial <- function(data, dists, fit, node, evidence, parents, predictions = NULL){
  if (dists[[node]] != "binomial"){
    stop("The node to predict should follow a binomial distribution.")
  }
  if (node %in% names(evidence)){
    # node is an evidence
    node_hat <- evidence[[node]]
  } else {
    if (length(parents)==0){
      # no parents
      node_hat <- predict_root(data, dists, node)
    } else {
      eq <- fit[[node]]
      eq_names <- colnames(eq) %||% names(eq)
      eq <- setNames(as.vector(eq), gsub(".*\\|", "", eq_names))

      bin.nodes <- intersect(names(dists)[which(dists %in% c("binomial","multinomial"))],parents)
      if (length(bin.nodes)>0){
        other.nodes <- parents[-which(parents %in% bin.nodes)]
      } else {
        other.nodes <- parents
      }

      predictions_tmp <- predictions[other.nodes]
      predictions_tmp <- lapply(predictions_tmp, function(l){
        l[1]
      })
      continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
      names(continuous_part) <- c()

      if (length(bin.nodes)>0){
        probabilities <- predictions[bin.nodes]

        bin.nodes.evidence <- intersect(names(evidence),bin.nodes)
        bin.nodes.constant <- names(which(lapply(predictions[bin.nodes],length)==1))
        if (length(bin.nodes.constant)>0){
          bin.nodes.evidence <- c(bin.nodes.evidence,bin.nodes.constant)
          bin.nodes.evidence <- unique(bin.nodes.evidence)
        }
        if (length(bin.nodes.evidence)>0){
          # at least one bin nodes is an evidence
          for (i in (1:length(bin.nodes.evidence))){
            proba_tmp <- rep(0,length(levels(probabilities[[bin.nodes.evidence[i]]])))
            names(proba_tmp) <- levels(data[[bin.nodes.evidence[i]]])
            proba_tmp[grep(probabilities[[bin.nodes.evidence[i]]],names(proba_tmp))] <- 1
            probabilities[[bin.nodes.evidence[i]]] <- proba_tmp
          }
        }
        names(probabilities) <- bin.nodes
        levels_list <- lapply(probabilities,function(x) names(x))
        combinations <- expand.grid(levels_list)

        bin_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "binomial")]
        multi_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "multinomial")]

        if(length(bin_parents) > 0) {
          mat_bin <- as.matrix(sapply(combinations[bin_parents], function(x) as.numeric(as.character(x))))
        } else {
          mat_bin <- NULL
        }

        if(length(multi_parents) > 0) {
          df_multi <- combinations[multi_parents]
          for(n in multi_parents) df_multi[[n]] <- factor(df_multi[[n]], levels = names(probabilities[[n]]))
          mat_multi <- model.matrix(~ ., data = df_multi)
          mat_multi <- mat_multi[, colnames(mat_multi) != "(Intercept)", drop = FALSE]
        } else {
          mat_multi <- NULL
        }

        dummy_matrix <- cbind(mat_bin, mat_multi)

        combinations_tmp <- dummy_matrix %*% eq[colnames(dummy_matrix)]

        proba_cond_values <- 1/(1+exp(- (continuous_part + rowSums(combinations_tmp))))

        prob_grid <- expand.grid(probabilities)
        combination_probabilities <- apply(prob_grid, 1, prod)

        node_hat <- sum(proba_cond_values * combination_probabilities)
      } else {
        node_hat <- 1/(1+exp(-continuous_part))
      }

      node_hat <- c(1 - node_hat, node_hat)
      names(node_hat) <- levels(data[[node]])
    }
  }
  return(node_hat)
}

#' Perform upstream inference with ABN
#'
#' Main function to predict the distribution of a multinomial node given its parents only
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param parents The parents of the node to predict.
#' @param predictions The estimated predictions of the upstream nodes  (can be empty if the parents of the node to predict are evidence).
#' @return The predicted distribution of the node of interest.
#' @import igraph
#' @export
#'
predict_node_from_parent_multinomial <- function(data, dists, fit, node, evidence, parents, predictions = NULL){
  if (dists[[node]] != "multinomial"){
    stop("The node to predict should follow a multinomial distribution.")
  }
  if (node %in% names(evidence)){
    # node is an evidence
    node_hat <- evidence[[node]]
  } else {
    if (length(parents)==0){
      # no parents
      node_hat <- predict_root(data, dists, node)
    } else {
      eq <- fit[[node]]
      eq_names <- colnames(eq) %||% names(eq)
      eq <- setNames(as.vector(eq), gsub(".*\\|", "", eq_names))

      bin.nodes <- intersect(names(dists)[which(dists %in% c("binomial","multinomial"))],parents)
      if (length(bin.nodes)>0){
        other.nodes <- parents[-which(parents %in% bin.nodes)]
      } else {
        other.nodes <- parents
      }

      predictions_tmp <- predictions[other.nodes]
      predictions_tmp <- lapply(predictions_tmp, function(l){
        l[1]
      })

      levels.multi <- levels(data[[node]])
      levels.cat <- levels.multi[2:length(levels.multi)]

      continuous_part <- lapply(levels.cat,function(level){
        eq[paste0("intercept.",level)] + sum(eq[paste0(other.nodes,".",level)]*unlist(predictions_tmp))
      })
      names(continuous_part) <- levels.cat

      if (length(bin.nodes)>0){
        probabilities <- predictions[bin.nodes]

        bin.nodes.evidence <- intersect(names(evidence),bin.nodes)
        if (length(bin.nodes.evidence)>0){
          # at least one bin nodes is an evidence
          for (i in (1:length(bin.nodes.evidence))){
            proba_tmp <- rep(0,length(levels(probabilities[[bin.nodes.evidence[i]]])))
            names(proba_tmp) <- levels(data[[bin.nodes.evidence[i]]])
            proba_tmp[grep(probabilities[[bin.nodes.evidence[i]]],names(proba_tmp))] <- 1
            probabilities[[bin.nodes.evidence[i]]] <- proba_tmp
          }
        }
        names(probabilities) <- bin.nodes
        levels_list <- lapply(probabilities,function(x) names(x))
        combinations <- expand.grid(levels_list)

        bin_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "binomial")]
        multi_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "multinomial")]

        if(length(bin_parents) > 0) {
          mat_bin <- as.matrix(sapply(combinations[bin_parents], function(x) as.numeric(as.character(x))))
        } else {
          mat_bin <- NULL
        }

        if(length(multi_parents) > 0) {
          df_multi <- combinations[multi_parents]
          for(n in multi_parents) df_multi[[n]] <- factor(df_multi[[n]], levels = names(probabilities[[n]]))
          mat_multi <- model.matrix(~ ., data = df_multi)
          mat_multi <- mat_multi[, colnames(mat_multi) != "(Intercept)", drop = FALSE]
        } else {
          mat_multi <- NULL
        }

        dummy_matrix <- cbind(mat_bin, mat_multi)

        numerators_matrix <- matrix(0, nrow = nrow(combinations), ncol = length(levels.cat))
        colnames(numerators_matrix) <- levels.cat
        for (level in levels.cat) {
          # Get discrete coefficients for this level
          target_names <- paste0(colnames(dummy_matrix), ".", level)
          eq_tmp <- eq[target_names]
          eq_tmp[is.na(eq_tmp)] <- 0

          # Linear predictor -> exp(eta)
          eta <- continuous_part[[level]] + as.vector(dummy_matrix %*% eq_tmp)
          numerators_matrix[, level] <- exp(eta)
        }
        scenario_denominators <- 1 + rowSums(numerators_matrix)
        proba_matrix <- numerators_matrix / scenario_denominators

        prob_grid <- expand.grid(probabilities)
        combination_probabilities <- apply(prob_grid, 1, prod)

        node_hat_others <- colSums(proba_matrix * combination_probabilities)
        node_hat_baseline <- sum((1 / scenario_denominators) * combination_probabilities)
        node_hat <- c(node_hat_baseline, node_hat_others)
        names(node_hat) <- levels.multi
      } else {
        numerator <- sapply(continuous_part, function(continuous_part_tmp){
          exp(continuous_part_tmp)
        })

        denominator <- 1+sum(numerator)

        node_hat_others <- numerator / denominator
        node_hat_baseline <- 1 / denominator

        node_hat <- c(node_hat_baseline, node_hat_others)
        names(node_hat) <- levels.multi
      }
    }
  }
  return(node_hat)
}

#' Perform downstream inference with ABN
#'
#' Main function to predict the distribution of a node given one of its Gaussian child and its parents
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param child A child of the node to predict.
#' @param parents The parents of the child.
#' @param predictions The estimated predictions of the downstream nodes (must contain at least a first prediction of the node to predict if the child and its parents are evidence).
#' @return The predicted distribution of the node of interest.
#' @import igraph
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#' graph <- igraph::graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#'
#' node <- "g1"
#' child <- find_children(graph, node)[2]
#' parents <- find_parents(graph, child)
#' evidence <- list("g2" = 3, "p1" = 3, "b2" = "y")
#' evidence <- check_evidence(data, mydists, hypothesis = node,evidence) # check the format of the evidence
#'
#' predictions <- list("g1" = c(0,1)) # a first estimate of the node g1
#'
#' predictions <- predict_node_from_children_gaussian(data, mydists, myfit, node, evidence, child, parents, predictions)
#' predictions  # the predicted distribution of g2
#' @export
#'
predict_node_from_children_gaussian <- function(data, dists, fit, node, evidence, child, parents, predictions){
  if (dists[[child]] != "gaussian"){
    stop("The child should follow a Gaussian distribution.")
  }
  if (!node %in% names(predictions)){
    stop("Predictions must contain at least a first prediction of the node to predict.")
  }
  if (node %in% names(evidence)) {
    return(predictions[[node]])
  }
  if (!all(c(child,parents) %in% names(predictions))){
    nodes <- setdiff(c(child,parents),names(predictions))
    if (!all(nodes %in% names(evidence))){
      stop("Not enough information about the downstream nodes.")
    }
    predictions <- c(predictions,evidence)
  }

  eq <- fit[[child]]
  eq_names <- colnames(eq) %||% names(eq)
  eq <- setNames(as.vector(eq), gsub(".*\\|", "", eq_names))

  bin.nodes <- intersect(names(dists)[which(dists %in% c("binomial","multinomial"))],parents)
  bin.nodes <- setdiff(bin.nodes,node)
  other.nodes <- setdiff(parents, bin.nodes)
  other.nodes <- setdiff(other.nodes, node)

  predictions_tmp <- predictions[other.nodes]
  predictions_tmp <- lapply(predictions_tmp, function(l){ l[1] })
  continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
  names(continuous_part) <- c()

  y_val <- predictions[[child]][1]

  if (length(predictions[[child]])>1){
    y_var <- predictions[[child]][2]
  } else {
    y_var = 1
  }

  compute_update <- function(intercept_tmp, node_type){
    switch(node_type,
           "binomial" = {
             p_prior <- predictions[[node]]
             if (length(p_prior)==1){
               ev_index <- as.numeric(as.factor(p_prior))
               p_vector <- c(0, 0)
               p_vector[ev_index] <- 1
               names(p_vector) <- levels(p_prior)
             } else {
               p_vector <- p_prior
             }
             numerator <- function(x) LogL_gaussian(y = y_val, x, coef = eq[[node]], var = y_var, intercept_tmp) + log(prior_binomial(x, p_vector[2]))

             log_num0 <- numerator(0)
             log_num1 <- numerator(1)
             max_log <- max(log_num0, log_num1)

             exp_num0 <- exp(log_num0 - max_log)
             exp_num1 <- exp(log_num1 - max_log)
             denominator <- exp_num0 + exp_num1
             prob_0 <- exp_num0 / denominator
             prob_1 <- 1 - prob_0
             results <- c(prob_0, prob_1)
             names(results) <- levels(data[[node]])
             return(results)
           },
           "gaussian" = {
             mu_prior <- predictions[[node]][1]
             sigma_prior <- predictions[[node]][2]

             build_integrand <- function(pow, shift = 0) {
               function(x) {
                 ((x - shift)^pow) *
                   L_gaussian(y = y_val, x, coef = eq[[node]], var = y_var, intercept_tmp) *
                   prior_gaussian(x, mu_prior, sigma_prior)
               }
             }

             run_integration <- function(f) {
               res <- try(integrate(f, -Inf, Inf)$value, silent = TRUE)
               if (inherits(res, "try-error")) {
                 low_b <- mu_prior - 10 * sigma_prior
                 upp_b <- mu_prior + 10 * sigma_prior
                 res <- try(integrate(f, low_b, upp_b)$value, silent = TRUE)
                 if (inherits(res, "try-error")) return(NA)
               }
               return(res)
             }

             denominator <- run_integration(build_integrand(pow = 0))
             if (denominator <= 0 || is.na(denominator)) {
               warning(paste0("Numerical underflow for Gaussian-Gaussian update at node ", node, ". Reverting to prior."))
               return(predictions[[node]])
             }
             f_mean <- build_integrand(1)
             m1 <- run_integration(f_mean) / denominator
             f_var <- build_integrand(2, shift = m1)
             posterior_variance <- run_integration(f_var) / denominator
             results <- c(m1, posterior_variance)
             return(results)
           },
           "poisson" = {
             lambda_prior <- predictions[[node]]
             max_x <- max(1000, 4 * lambda_prior)
             sum_val <- function(pow) {
               sum(sapply(0:max_x, function(x) (x^pow) * L_gaussian(y = y_val, x, coef = eq[[node]], var = y_var, intercept_tmp) *
                            prior_poisson(x, lambda_prior)))
             }
             denominator <- sum_val(0)
             if (denominator <= 0 || is.na(denominator)) return(lambda_prior)
             numerator <- sum_val(1)
             results <- numerator / denominator
             return(results)
           },
           "multinomial" = {
             p_prior <- predictions[[node]]
             if (length(p_prior) == 1) {
               levels_node <- levels(data[[node]])
               p_vector <- setNames(numeric(length(levels_node)), levels_node)
               p_vector[as.character(p_prior)] <- 1
             } else {
               p_vector <- p_prior
               levels_node <- names(p_prior)
             }

             unnormalized <- sapply(levels_node, function(l){
               coef_name <- paste0(node, l)
               node_coef <- if(coef_name %in% names(eq)) eq[[coef_name]] else 0
               lik <- L_gaussian(y = y_val, x = 1, coef = node_coef, var = y_var, intercept_tmp)
               return(lik * p_vector[l])
             })
             if (sum(unnormalized) == 0) return(p_vector)
             return(unnormalized / sum(unnormalized))
           }
    )
  }

  if (length(bin.nodes)>0){
    probabilities <- lapply(bin.nodes, function(bin) {
      lvls <- levels(data[[bin]])
      if (bin %in% names(evidence)) {
        vec <- as.numeric(lvls == evidence[[bin]])
        setNames(vec, lvls)
      } else {
        predictions[[bin]]
      }
    })
    names(probabilities) <- bin.nodes
    levels_list <- lapply(probabilities,function(x) names(x))
    combinations <- expand.grid(levels_list)

    bin_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "binomial")]
    multi_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "multinomial")]

    if(length(bin_parents) > 0) {
      mat_bin <- as.matrix(sapply(combinations[bin_parents], function(x) as.numeric(as.character(x))))
    } else {
      mat_bin <- NULL
    }

    if(length(multi_parents) > 0) {
      df_multi <- combinations[multi_parents]
      for(n in multi_parents) df_multi[[n]] <- factor(df_multi[[n]], levels = names(probabilities[[n]]))
      mat_multi <- model.matrix(~ ., data = df_multi)
      mat_multi <- mat_multi[, colnames(mat_multi) != "(Intercept)", drop = FALSE]
    } else {
      mat_multi <- NULL
    }

    dummy_matrix <- cbind(mat_bin, mat_multi)

    combinations_tmp <- dummy_matrix %*% eq[colnames(dummy_matrix)]

    proba_cond_values <- apply(combinations_tmp, 1, function(c) compute_update(continuous_part + c, dists[[node]]))

    prob_grid <- expand.grid(probabilities)
    comb_probs <- apply(prob_grid, 1, prod)

    if (is.matrix(proba_cond_values)) {
      final_res <- proba_cond_values %*% comb_probs
      final_res <- as.vector(final_res)
      if (dists[[node]] %in% c("binomial","multinomial")) names(final_res) <- levels(data[[node]])
    } else {
      final_res <- sum(proba_cond_values * comb_probs)
    }
  } else {
    final_res <- compute_update(continuous_part, dists[[node]])
  }
  return(final_res)
}

#' Perform downstream inference with ABN
#'
#' Main function to predict the distribution of a node given one of its Poisson child and its parents
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param child A child of the node to predict.
#' @param parents The parents of the child.
#' @param predictions The estimated predictions of the downstream nodes (must contain at least a first prediction of the node to predict if the child and its parents are evidence).
#' @return The predicted distribution of the node of interest.
#' @import igraph
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#' graph <- igraph::graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#'
#' node <- "b1"
#' child <- find_children(graph, node)[1]
#' parents <- find_parents(graph, child)
#' evidence <- list("p1" = 3, "p2" = 4)
#' evidence <- check_evidence(data, mydists, hypothesis = node,evidence) # check the format of the evidence
#'
#' predictions <- list("b1" = c(0.5,0.5)) # a first estimate of the node b1
#'
#' predictions <- predict_node_from_children_poisson(data, mydists, myfit, node, evidence, child, parents, predictions)
#' predictions  # the predicted distribution of g2
#' @export
#'
predict_node_from_children_poisson <- function(data, dists, fit, node, evidence, child, parents, predictions){
  if (dists[[child]] != "poisson"){
    stop("The child should follow a Poisson distribution.")
  }
  if (!node %in% names(predictions)){
    stop("Predictions must contain at least a first prediction of the node to predict.")
  }

  if (node %in% names(evidence)) {
    return(predictions[[node]])
  }
  if (!all(c(child,parents) %in% names(predictions))){
    nodes <- setdiff(c(child,parents),names(predictions))
    if (!all(nodes %in% names(evidence))){
      stop("Not enough information about the downstream nodes.")
    }
    predictions <- c(predictions,evidence)
  }

  eq <- fit[[child]]
  eq_names <- colnames(eq) %||% names(eq)
  eq <- setNames(as.vector(eq), gsub(".*\\|", "", eq_names))

  bin.nodes <- intersect(names(dists)[which(dists %in% c("binomial","multinomial"))],parents)
  bin.nodes <- setdiff(bin.nodes,node)
  other.nodes <- setdiff(parents, bin.nodes)
  other.nodes <- setdiff(other.nodes, node)

  predictions_tmp <- predictions[other.nodes]
  predictions_tmp <- lapply(predictions_tmp,function(l){
    l[1]
  })
  continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
  names(continuous_part) <- c()

  y_val <- predictions[[child]]

  compute_update <- function(intercept_tmp, node_type){
    switch(node_type,
           "binomial" = {
             p_prior <- predictions[[node]]
             if (length(p_prior)==1){
               ev_index <- as.numeric(as.factor(p_prior))
               p_vector <- c(0, 0)
               p_vector[ev_index] <- 1
               names(p_vector) <- levels(p_prior)
             } else {
               p_vector <- p_prior
             }
             numerator <- function(x) exp(LogL_poisson(y = y_val, x, coef = eq[[node]],intercept_tmp)) * prior_binomial(x, p_vector[2])
             denominator <- numerator(0) + numerator(1)

             if (is.na(denominator) || denominator == 0) return(p_vector) # Numerical fallback
             results <- c(numerator(0)/denominator, 1-numerator(0)/denominator)
             names(results) <- levels(data[[node]])
             return(results)
           },
           "gaussian" = {
             mu_prior <- predictions[[node]][1]
             sigma_prior <- predictions[[node]][2]

             build_integrand <- function(pow, shift = 0) {
               function(x) {
                 ((x - shift)^pow) *
                   exp(LogL_poisson(y = y_val, x, coef = eq[[node]], intercept_tmp)) *
                   prior_gaussian(x, mu_prior, sigma_prior)
               }
             }

             run_integration <- function(f) {
               res <- try(integrate(f, -Inf, Inf)$value, silent = TRUE)
               if (inherits(res, "try-error")) {
                 res <- integrate(f, mu_prior - 5 * sigma_prior, mu_prior + 5 * sigma_prior, rel.tol=1e-6)$value
               }
               return(res)
             }

             denominator <- run_integration(build_integrand(0))
             if (denominator <= 0 || is.na(denominator)) {
               warning(paste0("Numerical underflow for node ", node, ". Reverting to prior."))
               return(predictions[[node]])
             }
             f_mean <- build_integrand(1)
             m1 <- run_integration(f_mean) / denominator
             f_var <- build_integrand(2, shift = m1)
             posterior_variance <- run_integration(f_var) / denominator

             results <- c(m1, posterior_variance)

             if (is.na(posterior_variance) || posterior_variance <= 0) {
               warning(paste0("Numerical issues with variance of node ", node, ". Reverting to prior."))
               return(predictions[[node]])
             }
             return(results)
           },
           "poisson" = {
             lambda_prior <- predictions[[node]]
             max_x <- max(1000, 4 * lambda_prior)
             x_vals <- 0:max_x
             log_posterior_vec <- sapply(x_vals, function(x) {
               log(prior_poisson(x, lambda_prior)) + LogL_poisson(y = y_val, x, coef = eq[[node]], intercept_tmp)
             })
             log_posterior_vec[is.na(log_posterior_vec)] <- -Inf
             M <- max(log_posterior_vec)
             if (is.infinite(M)) return(lambda_prior)
             weights <- exp(log_posterior_vec - M)
             denominator <- sum(weights)
             numerator  <- sum(x_vals * weights)
             results <- numerator / denominator
             return(results)
           },
           "multinomial" = {
             p_prior <- predictions[[node]]
             if (length(p_prior) == 1) {
               levels_node <- levels(data[[node]])
               p_vector <- setNames(numeric(length(levels_node)), levels_node)
               p_vector[as.character(p_prior)] <- 1
             } else {
               p_vector <- p_prior
               levels_node <- names(p_prior)
             }

             unnormalized <- sapply(levels_node, function(l){
               coef_name <- paste0(node, l)
               node_coef <- if(coef_name %in% names(eq)) eq[[coef_name]] else 0 # to verify in practice
               log_lik <- LogL_poisson(y = y_val, x = 1, coef = node_coef, intercept_tmp)
               return(exp(log_lik) * p_vector[l])
             })
             if (sum(unnormalized) == 0) return(p_vector)
             return(unnormalized / sum(unnormalized))
           }
    )
  }

  if (length(bin.nodes)>0){
    probabilities <- lapply(bin.nodes, function(bin) {
      lvls <- levels(data[[bin]])
      if (bin %in% names(evidence)) {
        vec <- as.numeric(lvls == evidence[[bin]])
        setNames(vec, lvls)
      } else {
        predictions[[bin]]
      }
    })
    names(probabilities) <- bin.nodes
    levels_list <- lapply(probabilities,function(x) names(x))
    combinations <- expand.grid(levels_list)

    bin_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "binomial")]
    multi_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "multinomial")]

    if(length(bin_parents) > 0) {
      mat_bin <- as.matrix(sapply(combinations[bin_parents], function(x) as.numeric(as.character(x))))
    } else {
      mat_bin <- NULL
    }

    if(length(multi_parents) > 0) {
      df_multi <- combinations[multi_parents]
      for(n in multi_parents) df_multi[[n]] <- factor(df_multi[[n]], levels = names(probabilities[[n]]))
      mat_multi <- model.matrix(~ ., data = df_multi)
      mat_multi <- mat_multi[, colnames(mat_multi) != "(Intercept)", drop = FALSE]
    } else {
      mat_multi <- NULL
    }

    dummy_matrix <- cbind(mat_bin, mat_multi)

    combinations_tmp <- dummy_matrix %*% eq[colnames(dummy_matrix)]

    proba_cond_values <- apply(combinations_tmp, 1, function(c) compute_update(continuous_part + c, dists[[node]]))

    prob_grid <- expand.grid(probabilities)
    comb_probs <- apply(prob_grid, 1, prod)

    if (is.matrix(proba_cond_values)) {
      final_res <- proba_cond_values %*% comb_probs
      final_res <- as.vector(final_res)
      if (dists[[node]] %in% c("binomial","multinomial")) names(final_res) <- levels(data[[node]])
    } else {
      final_res <- sum(proba_cond_values * comb_probs)
    }
  } else {
    final_res <- compute_update(continuous_part, dists[[node]])
  }
  return(final_res)
}

#' Perform downstream inference with ABN
#'
#' Main function to predict the distribution of a node given one of its binomial child and its parents
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param child A child of the node to predict.
#' @param parents The parents of the child.
#' @param predictions The estimated predictions of the downstream nodes (must contain at least a first prediction of the node to predict if the child and its parents are evidence).
#' @return The predicted distribution of the node of interest.
#' @import igraph
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#' graph <- igraph::graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#'
#' node <- "g1"
#' child <- find_children(graph, node)[1]
#' parents <- find_parents(graph, child)
#' evidence <- list("b3" = "y",  "b1" = "y", "b2" = "y")
#' evidence <- check_evidence(data, mydists, hypothesis = node,evidence) # check the format of the evidence
#'
#' predictions <- list("g1" = c(0,1)) # a first estimate of the node g1
#'
#' predictions <- predict_node_from_children_binomial(data, mydists, myfit, node, evidence, child, parents, predictions)
#' predictions  # the predicted distribution of g2
#' @export
#'
predict_node_from_children_binomial <- function(data, dists, fit, node, evidence, child, parents, predictions){
  if (dists[[child]] != "binomial"){
    stop("The child should follow a binomial distribution.")
  }
  if (!node %in% names(predictions)){
    stop("Predictions must contain at least a first prediction of the node to predict.")
  }
  if (node %in% names(evidence)) {
    return(predictions[[node]])
  }
  if (!all(c(child,parents) %in% names(predictions))){
    nodes <- setdiff(c(child,parents),names(predictions))
    if (!all(nodes %in% names(evidence))){
      stop("Not enough information about the downstream nodes.")
    }
    predictions <- c(predictions,evidence)
  }

  eq <- fit[[child]]
  eq_names <- colnames(eq) %||% names(eq)
  eq <- setNames(as.vector(eq), gsub(".*\\|", "", eq_names))

  bin.nodes <- intersect(names(dists)[which(dists %in% c("binomial","multinomial"))],parents)
  bin.nodes <- setdiff(bin.nodes,node)
  other.nodes <- setdiff(parents, bin.nodes)
  other.nodes <- setdiff(other.nodes, node)

  predictions_tmp <- predictions[other.nodes]
  predictions_tmp <- lapply(predictions_tmp, function(l){ l[1] })
  continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
  names(continuous_part) <- c()

  if (child %in% names(evidence)){
    p_pred <- c(0,0)
    names(p_pred) <- levels(data[[child]])
    p_pred[names(p_pred) == evidence[[child]]] <- 1
  } else {
    p_pred <- predictions[[child]]
  }

  compute_update <- function(intercept_tmp, node_type){
    switch(node_type,
        "binomial" = {
           p_prior <- predictions[[node]]
           if (length(p_prior)==1){
             ev_index <- as.numeric(as.factor(p_prior))
             p_vector <- c(0, 0)
             p_vector[ev_index] <- 1
             names(p_vector) <- levels(p_prior)
           } else {
             p_vector <- p_prior
           }
           numerator <- function(x, y) L_binomial(y, x, coef = eq[[node]], intercept_tmp) * prior_binomial(x, p_vector[2])
           denominator <- (numerator(0, 0) + numerator(1, 0)) * p_pred[1] + (numerator(0, 1) + numerator(1, 1)) * p_pred[2]
           prob_0 <- (numerator(0, 0) * p_pred[1] + numerator(0, 1) * p_pred[2]) / denominator
           results <- c(prob_0, 1 - prob_0)
           names(results) <- levels(data[[node]])
           return(results)
        },
        "gaussian" = {
          mu_prior <- predictions[[node]][1]
          sigma_prior <- predictions[[node]][2]

          build_integrand <- function(pow, shift = 0) {
            function(x) {
              total_lik_density <- L_binomial(y = 0, x, coef = eq[[node]], intercept_tmp) * p_pred[1] +
                L_binomial(y = 1, x, coef = eq[[node]], intercept_tmp) * p_pred[2]
              return(((x - shift)^pow) * total_lik_density * prior_gaussian(x, mu_prior, sigma_prior))
            }
          }

          run_integration <- function(f) {
            res <- try(integrate(f, -Inf, Inf)$value, silent = TRUE)
            if (inherits(res, "try-error")) {
              low_b <- mu_prior - 10 * sigma_prior
              upp_b <- mu_prior + 10 * sigma_prior
              res <- try(integrate(f, low_b, upp_b)$value, silent = TRUE)
              if (inherits(res, "try-error")) return(NA)
            }
            return(res)
          }
          denominator <- run_integration(build_integrand(pow = 0))

          if (denominator <= 0 || is.na(denominator)) {
            warning(paste0("Numerical underflow for Gaussian-Binomial update at node ", node, ". Reverting to prior."))
            return(predictions[[node]])
          }
          m1 <- run_integration(build_integrand(pow = 1)) / denominator
          posterior_variance <- run_integration(build_integrand(pow = 2, shift = m1)) / denominator

          if (is.na(posterior_variance) || posterior_variance <= 0) {
            warning(paste0("Numerical issues with variance calculation at node ", node, ". Reverting to prior."))
            return(predictions[[node]])
          }
          results <- c(m1, posterior_variance)
          return(results)
        },
        "poisson" = {
          lambda_prior <- predictions[[node]]
          max_x <- max(1000, 4 * lambda_prior)
          sum_val <- function(pow, y) {
            sum(sapply(0:max_x, function(x) (x^pow) * L_binomial(y, x, coef = eq[[node]], intercept_tmp) *
                         prior_poisson(x, lambda_prior) * p_pred[y + 1]))
          }
          denominator <- sum_val(0, 0) + sum_val(0, 1)
          if (denominator <= 0 || is.na(denominator)) return(lambda_prior)
          numerator <- sum_val(1, 0) + sum_val(1, 1)
          results <- numerator / denominator
          return(results)
        },
        "multinomial" = {
          p_prior <- predictions[[node]]
          if (length(p_prior) == 1) {
            levels_node <- levels(data[[node]])
            p_vector <- setNames(numeric(length(levels_node)), levels_node)
            p_vector[as.character(p_prior)] <- 1
          } else {
            p_vector <- p_prior
            levels_node <- names(p_prior)
          }
          unnormalized <- sapply(levels_node, function(l){
            coef_name <- paste0(node, l)
            node_coef <- if(coef_name %in% names(eq)) as.numeric(eq[[coef_name]]) else 0
            lik <- L_binomial(y = 0, x = 1, coef = node_coef, intercept_tmp)*p_pred[1]+L_binomial(y = 1, x = 1, coef = node_coef, intercept_tmp)*p_pred[2]
            return(lik * p_vector[l])
          })

          if (sum(unnormalized) == 0) return(p_vector)
          return(unnormalized / sum(unnormalized))
        }
      )
  }

  if (length(bin.nodes)>0){
    probabilities <- lapply(bin.nodes, function(bin) {
      lvls <- levels(data[[bin]])
      if (bin %in% names(evidence)) {
        setNames(as.numeric(lvls == evidence[[bin]]), lvls)
      } else {
        predictions[[bin]]
      }
    })
    names(probabilities) <- bin.nodes

    levels_list <- lapply(probabilities,function(x) names(x))
    combinations <- expand.grid(levels_list)

    bin_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "binomial")]
    multi_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "multinomial")]

    if(length(bin_parents) > 0) {
      mat_bin <- as.matrix(sapply(combinations[bin_parents], function(x) as.numeric(as.character(x))))
    } else {
      mat_bin <- NULL
    }

    if(length(multi_parents) > 0) {
      df_multi <- combinations[multi_parents]
      for(n in multi_parents) df_multi[[n]] <- factor(df_multi[[n]], levels = names(probabilities[[n]]))
      mat_multi <- model.matrix(~ ., data = df_multi)
      mat_multi <- mat_multi[, colnames(mat_multi) != "(Intercept)", drop = FALSE]
    } else {
      mat_multi <- NULL
    }

    dummy_matrix <- cbind(mat_bin, mat_multi)

    combinations_tmp <- dummy_matrix %*% eq[colnames(dummy_matrix)]

    proba_cond_values <- apply(combinations_tmp, 1, function(c) compute_update(continuous_part + c, dists[[node]]))

    prob_grid <- expand.grid(probabilities)
    comb_probs <- apply(prob_grid, 1, prod)

    if (is.matrix(proba_cond_values)) {
      final_res <- proba_cond_values %*% comb_probs
      final_res <- as.vector(final_res)
      if (dists[[node]] %in% c("binomial","multinomial")) names(final_res) <- levels(data[[node]])
    } else {
      final_res <- sum(proba_cond_values * comb_probs)
    }
  } else {
    final_res <- compute_update(continuous_part, dists[[node]])
  }
  return(final_res)
}

#' Perform downstream inference with ABN
#'
#' Main function to predict the distribution of a node given one of its multinomial child and its parents
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param fit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict.
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param child A child of the node to predict.
#' @param parents The parents of the child.
#' @param predictions The estimated predictions of the downstream nodes (must contain at least a first prediction of the node to predict if the child and its parents are evidence).
#' @return The predicted distribution of the node of interest.
#' @import igraph
#' @export
#'
predict_node_from_children_multinomial <- function(data, dists, fit, node, evidence, child, parents, predictions){
  if (dists[[child]] != "multinomial"){
    stop("The child should follow a multinomial distribution.")
  }
  if (!node %in% names(predictions)){
    stop("Predictions must contain at least a first prediction of the node to predict.")
  }
  if (node %in% names(evidence)) {
    return(predictions[[node]])
  }
  if (!all(c(child,parents) %in% names(predictions))){
    nodes <- setdiff(c(child,parents),names(predictions))
    if (!all(nodes %in% names(evidence))){
      stop("Not enough information about the downstream nodes.")
    }
    predictions <- c(predictions,evidence)
  }

  eq <- fit[[child]]
  eq_names <- colnames(eq) %||% names(eq)
  eq <- setNames(as.vector(eq), gsub(".*\\|", "", eq_names))

  bin.nodes <- intersect(names(dists)[which(dists %in% c("binomial","multinomial"))],parents)
  bin.nodes <- setdiff(bin.nodes,node)
  other.nodes <- setdiff(parents, bin.nodes)
  other.nodes <- setdiff(other.nodes, node)

  predictions_tmp <- predictions[other.nodes]
  predictions_tmp <- lapply(predictions_tmp,function(l){
    l[1]
  })

  lvl_names <- levels(data[[child]])

  continuous_part <- sapply(lvl_names, function(p) {
    int_name <- paste0("intercept.", p)

    if (int_name %in% names(eq)) {
      val <- eq[int_name] + sum(eq[paste0(other.nodes, ".", p)] * unlist(predictions_tmp))
      return(val)
    } else {
      return(0)
    }
  })
  names(continuous_part) <- lvl_names

  if (child %in% names(evidence)){
    p_pred <- rep(0,length(lvl_names))
    names(p_pred) <- lvl_names
    p_pred[names(p_pred) == evidence[[child]]] <- 1
  } else {
    p_pred <- predictions[[child]]
  }

  node_coef_vec <- sapply(lvl_names, function(lvl) {
    c_name <- paste0(node, ".", lvl)
    if(c_name %in% names(eq)) as.numeric(eq[[c_name]]) else 0
  })

  compute_update <- function(intercept_tmp, node_type){
    switch(node_type,
           "binomial" = {
             p_prior <- predictions[[node]]

             if (length(p_prior)==1){
               ev_index <- as.numeric(as.factor(p_prior))
               p_vector <- c(0, 0)
               p_vector[ev_index] <- 1
               names(p_vector) <- levels(p_prior)
             } else {
               p_vector <- p_prior
             }

             numerator <- function(x, y){
               L_multinomial(y, x, coef = node_coef_vec, intercept_tmp) * prior_binomial(x, p_vector[2])
             }
             denominator <- sum(sapply(seq_along(lvl_names), function(p){
               (numerator(0, p) +numerator(1, p)) * p_pred[p]
             }))
             prob_0 <- sum(sapply(seq_along(lvl_names), function(p){
               numerator(0, p) * p_pred[p]
             })) / denominator
             results <- c(prob_0, 1 - prob_0)
             names(results) <- levels(data[[node]])
             return(results)
           },
           "gaussian" = {
             mu_prior <- predictions[[node]][1]
             sigma_prior <- predictions[[node]][2]

             build_integrand <- function(pow, shift = 0) {
               function(x) {
                 total_lik_density <- numeric(length(x))

                 for (p in seq_along(lvl_names)) {
                   total_lik_density <- total_lik_density +
                     L_multinomial(y = p, x, coef = node_coef_vec, intercept_tmp) * p_pred[p]
                 }

                 # Combined density at coordinate points
                 return(((x - shift)^pow) * total_lik_density * prior_gaussian(x, mu_prior, sigma_prior))
               }
             }

             run_integration <- function(f) {
               res <- try(integrate(f, -Inf, Inf)$value, silent = TRUE)
               if (inherits(res, "try-error")) {
                 low_b <- mu_prior - 10 * sigma_prior
                 upp_b <- mu_prior + 10 * sigma_prior
                 res <- integrate(f, low_b, upp_b)$value
               }
               return(res)
             }

             denominator <- run_integration(build_integrand(pow = 0))
             if (denominator <= 0 || is.na(denominator)) {
               warning(paste0("Numerical underflow for Gaussian-Multinomial update at node ", node, ". Reverting to prior."))
               return(predictions[[node]])
             }
             m1 <- run_integration(build_integrand(pow = 1)) / denominator
             posterior_variance <- run_integration(build_integrand(pow = 2, shift = m1)) / denominator
             if (is.na(posterior_variance) || posterior_variance <= 0) {
               warning(paste0("Numerical issues with variance calculation at node ", node, ". Reverting to prior."))
               return(predictions[[node]])
             }
             results <- c(m1, posterior_variance)
             return(results)
           },
           "poisson" = {
             lambda_prior <- predictions[[node]]
             max_x <- max(1000, 4 * lambda_prior)

             sum_val <- function(pow, y) {
               sum(sapply(0:max_x, function(x) (x^pow) * L_binomial(y, x, coef = node_coef_vec, intercept_tmp) *
                            prior_poisson(x, lambda_prior) * p_pred[y]))
             }
             denominator <- sum(sapply(seq_along(lvl_names), function(p){
               sum_val(0, p)
             }))
             numerator <- sum(sapply(seq_along(lvl_names), function(p){ sum_val(1, p) }))
             results <- numerator / denominator
             return(results)
           },
           "multinomial" = {
             p_prior <- predictions[[node]]
             if (length(p_prior) == 1) {
               levels_node <- levels(data[[node]])
               p_vector <- setNames(numeric(length(levels_node)), levels_node)
               p_vector[as.character(p_prior)] <- 1
             } else {
               p_vector <- p_prior
               levels_node <- names(p_prior)
             }
             unnormalized <- sapply(levels_node, function(l){
               current_parent_coefs <- sapply(lvl_names, function(clvl) {
                 c_name <- paste0(node, l, ".", clvl) # to check in practice
                 if(c_name %in% names(eq)) as.numeric(eq[[c_name]]) else 0
               })

               lik <- sum(sapply(seq_along(lvl_names),function(p){
                 L_binomial(y = p, x = 1, coef = current_parent_coefs, intercept_tmp) * p_pred[p]
               }))
               return(lik * p_vector[l])
             })

             if (sum(unnormalized) == 0) return(p_vector)
             return(unnormalized / sum(unnormalized))
           }
    )
  }

  if (length(bin.nodes)>0){
    probabilities <- lapply(bin.nodes, function(bin) {
      lvls <- levels(data[[bin]])
      if (bin %in% names(evidence)) {
        vec <- as.numeric(lvls == evidence[[bin]])
        setNames(vec, lvls)
      } else {
        predictions[[bin]]
      }
    })
    names(probabilities) <- bin.nodes
    levels_list <- lapply(probabilities,function(x) names(x))
    combinations <- expand.grid(levels_list)

    bin_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "binomial")]
    multi_parents <- names(combinations)[sapply(names(combinations), function(x) dists[[x]] == "multinomial")]

    if(length(bin_parents) > 0) {
      mat_bin <- as.matrix(sapply(combinations[bin_parents], function(x) as.numeric(as.character(x))))
    } else {
      mat_bin <- NULL
    }

    if(length(multi_parents) > 0) {
      df_multi <- combinations[multi_parents]
      for(n in multi_parents) df_multi[[n]] <- factor(df_multi[[n]], levels = names(probabilities[[n]]))
      mat_multi <- model.matrix(~ ., data = df_multi)
      mat_multi <- mat_multi[, colnames(mat_multi) != "(Intercept)", drop = FALSE]
    } else {
      mat_multi <- NULL
    }

    dummy_matrix <- cbind(mat_bin, mat_multi)

    scenario_effects <- sapply(lvl_names, function(lvl) {
      target_names <- paste0(colnames(dummy_matrix), ".", lvl)
      exists_mask <- target_names %in% names(eq)

      if (!any(exists_mask)) return(numeric(nrow(dummy_matrix)))

      relevant_eq <- as.numeric(eq[target_names[exists_mask]])
      relevant_mat <- dummy_matrix[, colnames(dummy_matrix)[exists_mask], drop = FALSE]

      as.vector(relevant_mat %*% relevant_eq)
    })
    colnames(scenario_effects) <- lvl_names

    proba_cond_values <- apply(scenario_effects, 1, function(c) {
      compute_update(continuous_part + c, dists[[node]])
    })
    prob_grid <- expand.grid(probabilities)
    comb_probs <- apply(prob_grid, 1, prod)

    if (is.matrix(proba_cond_values)) {
      final_res <- proba_cond_values %*% comb_probs
      final_res <- as.vector(final_res)
      if (dists[[node]] %in% c("binomial","multinomial")) names(final_res) <- levels(data[[node]])
    } else {
      final_res <- sum(proba_cond_values * comb_probs)
    }
  } else {
    final_res <- compute_update(continuous_part, dists[[node]])
  }
  return(final_res)
}

#' Find the parents of a node in a graph
#'
#' Find the parents of a given node in a given graph.
#'
#' @param graph A directed graph (igraph object).
#' @param node A node of the graph (either a node label or a number)
#' @return A vector containing the parents of the node (labels of the nodes if the graph is labeled)
#' @import igraph
#' @examples
#' g <- igraph::make_graph("Zachary") # undirected graph
#' find_parents(g, node = 1) # neighbors of the node 1
#' @export
find_parents <- function(graph,node){
  if (is_directed(graph)==FALSE){
    warning(paste0("The provided graph is not directed, this function will output the neighbors of the node ",node,"."))
  }
  if (!is_igraph(graph)){
    stop("The provided graph should be an igraph object.")
  }
  if (is.null(V(graph)$name)){
    if (is.character(node)){
      stop("The graph is not labeled, provide a number corresponding to one node.")
    } else if (node > length(V(graph))){
      stop(paste0("There are less than ", length(V(g)), " nodes in the graph."))
    }
    parents <- as.character(neighbors(graph, v= node,mode="in"))
  } else {
    if (is.character(node)){
      if (!node %in% V(graph)$name){
        stop("The node you selected is not part of the graph.")
      }
    } else {
      if (node > length(V(graph))){
        stop(paste0("There are less than ", length(V(g)), " nodes in the graph."))
      }
    }
    parents <- names(neighbors(graph, v= node,mode="in"))
  }
  return(parents)
}

#' Find the children of a node in a graph
#'
#' Find the children of a given node in a given graph.
#'
#' @param graph A directed graph (igraph object).
#' @param node A node of the graph (either a node label or a number)
#' @return A vector containing the children of the node (labels of the nodes if the graph is labeled)
#' @import igraph
#' @examples
#' g <- igraph::make_graph("Zachary") # undirected graph
#' find_children(g, node = 1) # neighbors of the node 1
#' @export
find_children <- function(graph,node){
  if (is_directed(graph)==FALSE){
    warning(paste0("The provided graph is not directed, this function will output the neighbors of the node ",node,"."))
  }
  if (!is_igraph(graph)){
    stop("The provided graph should be an igraph object.")
  }
  if (is.null(V(graph)$name)){
    if (is.character(node)){
      stop("The graph is not labeled, provide a number corresponding to one node.")
    } else if (node > length(V(graph))){
      stop(paste0("There are less than ", length(V(g)), " nodes in the graph."))
    }
    children <- as.character(neighbors(graph, v= node,mode="out"))
  } else {
    if (is.character(node)){
      if (!node %in% V(graph)$name){
        stop("The node you selected is not part of the graph.")
      }
    } else {
      if (node > length(V(graph))){
        stop(paste0("There are less than ", length(V(g)), " nodes in the graph."))
      }
    }
    children <- names(neighbors(graph, v= node,mode="out"))
  }
  return(children)
}

#' Predict a root node
#'
#' Predict a root node in the graph
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param node Root node to predict.
#' @return The predicted distribution of the root node.
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' prediction <- predict_root(data, mydists, node = "g3")
#' @export
predict_root <- function(data, dists, node) {
  switch(dists[[node]],
         "gaussian"  = mean(data[[node]]),
         "poisson"   = mean(data[[node]]),
         "binomial"  = prop.table(table(data[[node]])),
         "multinomial" = prop.table(table(data[[node]])),
         stop("Unsupported distribution type for root node ", dists[[node]])
  )
}

#' Compute the log-likelihood of a Poisson variable
#'
#' This function computes the log-likelihood of a Poisson-distributed variable
#'
#' @param y The observed count for the Poisson-distributed variable (must be non-negative).
#' @param x A value for a parent of the Poisson variable.
#' @param coef The coefficient that links `y` to `x`.
#' @param continuous_part The intercept and any additional contributions from other parents of `y`.
#' @return The value of the log-likelihood.
#' @export
LogL_poisson <- function(y, x, coef, continuous_part){
  lambda <- exp(coef*x + continuous_part)
  #dpois(y, lambda)
  #(lambda^y * exp(-lambda)) / gamma(y + 1)
  y * log(lambda) - lambda - lgamma(y + 1)
}

#' Compute Gaussian prior density
#'
#' This function calculates the density of a Gaussian distribution
#'
#' @param x The value at which the density is evaluated.
#' @param mu The expected value of the Gaussian distribution.
#' @param sigma2 The variance of the Gaussian distribution.
#' @return The density of the normal distribution at x.
#' @export
prior_gaussian <- function(x, mu, sigma2){
  dnorm(x, mean=mu, sd=sqrt(sigma2))
}

#' Compute the likelihood of a Gaussian variable
#'
#' This function computes the likelihood of a Gaussian-distributed variable
#'
#' @param y The observed value for the Gaussian-distributed variable.
#' @param x A value for a parent of the Gaussian variable.
#' @param coef The coefficient that links `y` to `x`.
#' @param var The variance of the Gaussian variable.
#' @param continuous_part The rest of the equation (intercept and links between the parents of the child and the child).
#' @return The value of the likelihood.
#' @export
L_gaussian <- function(y, x, coef, var, continuous_part){
  mu <-  coef*x + continuous_part
  sigma2 <- var
  dnorm(y, mean = mu, sd=sqrt(sigma2))
}

LogL_gaussian <- function(y, x, coef, var, continuous_part) {
  mu <- continuous_part + coef * x
  # Setting log = TRUE calculates the log-exponent directly without hitting the underflow wall!
  return(dnorm(y, mean = mu, sd = sqrt(var), log = TRUE))
}

#' Compute prior probability for a binomial prior distribution
#'
#' This function calculates the probability mass function of a binomial distribution
#'
#' @param x The observed outcome (0 or 1).
#' @param p The probability of sucess (between 0 and 1).
#' @return The probability of observing x given sucess probability p.
#' @export
prior_binomial <- function(x, p){
  dbinom(x, size=1, prob=p)
}

#' Compute the likelihood of a binomial variable
#'
#' This function computes the likelihood of a binomial-distributed variable
#'
#' @param y The observed binary outcome (numeric 0 or 1)
#' @param x A value for a parent of the binomial variable.
#' @param coef The coefficient that links `y` to `x`.
#' @param continuous_part The rest of the equation (intercept and links between the parents of the child and the child).
#' @return The value of the likelihood.
#' @export
L_binomial <- function(y, x, coef, continuous_part){
  mu <- plogis(continuous_part + coef * x)
  if (y == 1) {
    return(mu)
  } else {
    return(1 - mu)
  }
}

#' Compute the likelihood of a multinomial variable
#'
#' This function computes the likelihood of a multinomial-distributed variable
#'
#' @param y The observed multivariate outcome
#' @param x A value for a parent of the multinomial variable.
#' @param coef The coefficient that links `y` to `x`.
#' @param continuous_part The rest of the equation (intercept and links between the parents of the child and the child).
#' @return The value of the likelihood.
#' @export
L_multinomial <- function(y, x, coef, continuous_part){
  numerator <- function(y){
    exp(continuous_part[y] + coef[y] *x)
  }

  denominator <- 1 + sum(sapply(2:length(coef), function(p){
    numerator(p)}))

  if (y==1){
    return(1/denominator)
  } else {
    return(numerator(y)/denominator)
  }
}

#' Compute Prior Probability for a Poisson Distribution
#'
#' This function calculates the probability mass function of a Poisson distribution
#'
#' @param x The observed count (number of events, must be non-negative).
#' @param lambda The expected number of events (rate parameter, must be positive).
#' @return The probability of observing `x` events given rate `lambda`.
#' @export
prior_poisson <- function(x,lambda){
  dpois(x,lambda)
}

#' Plot ABN fitted network
#'
#' This function plots the ABN network with an emphasis on a specific node
#'
#' @param dag An adjacency matrix (can be the output of the function mostProbable()).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param node The node of interest.
#' @param order up or down indicating whether we predict the node from upstream or downstream.
#' @return A graph object.
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#'
#' node <- "b4"
#' order <- "up"
#' plot_node <- plot_Abn(dag, mydists, node, order)
#' plot_node
#' @import igraph
#' @export
plot_Abn <- function (dag, dists, node, order){
  if (!all(colnames(dag) %in% names(dists))){
    stop("The names of the nodes in dag do not correspond to the ones in dists.")
  }
  if (!node %in% colnames(dag)){
    stop("Choose a node that belongs to the network.")
  }
  if (! order %in% c("up","down")){
    stop("Choose either up or down as an argument for order.")
  }
  dists <- dists[colnames(dag)]
  name <- names(dists)

  graph <- graph_from_adjacency_matrix(t(dag))

  am.graph <- new(Class = "graphAM", adjMat = dag, edgemode = "directed")

  node.shape <- rep(c("circle", "box","ellipse","diamond"), 4)
  shape <- rep(node.shape[1], length(dists))
  shape[dists == "binomial"] <- node.shape[2]
  shape[dists == "poisson"] <- node.shape[3]
  shape[dists == "multinomial"] <- node.shape[4]
  names(shape) <- names(dists)

  node.fillcolor = c("lightblue", "brown3", "chartreuse3","chartreuse4")
  fillcolor <- rep(node.fillcolor[1], length(dists))
  names(fillcolor) <- names(dists)
  if (!is.null(node)) {
    markov.blanket <- abn::mb(dag, node = node,
                         data.dists = dists)
    fillcolor[names(dists) %in% node] <- node.fillcolor[2] # node in red

    if (order == "up"){
      parents <- find_parents(graph, node)
      fillcolor[names(dists) %in% parents] <- node.fillcolor[3]
    } else if (order == "down"){
      children <- find_children(graph, node)
      parents <- unique(unlist(sapply(children,function(l){find_parents(graph, l)})))
      fillcolor[names(dists) %in% setdiff(parents,node)] <- node.fillcolor[4]
      fillcolor[names(dists) %in% children] <- node.fillcolor[3]
    }
  }

  names.edges <- names(Rgraphviz::buildEdgeList(am.graph))
  edge.label <- rep(" ", length(names.edges))
  names(edge.label) <- names.edges
  edge.lwd <- rep(1, length(names.edges))
  class(edge.lwd) <- "character"
  names(edge.lwd) <- names.edges

  attrs <- list(graph = list(rankdir = "BT"), node = list(fontsize = 12,
                                                          fixedsize = FALSE), edge = list(arrowsize = 0.6,
                                                                                          color = "black", lty = "solid", fontsize = 12))
  nodeAttrs <- list(fillcolor = fillcolor, shape = shape)
  edgeAttrs <- list(label = edge.label, lwd = edge.lwd)
  am.graph <- layoutGraph(am.graph, attrs = attrs, nodeAttrs = nodeAttrs,
                          edgeAttrs = edgeAttrs)

  edgeRenderInfo(am.graph) <- list(arrowtail = "open")
  edgeRenderInfo(am.graph) <- list(arrowhead = "none")
  edgeRenderInfo(am.graph) <- list(lwd = edge.lwd)
  return(graph=am.graph)
}

#' Plot the procedure's workflow
#'
#' This function plots and saves an animated gif that represents the procedure's workflow.
#'
#' @param dag An adjacency matrix (can be the output of the function mostProbable()).
#' @param dists A list containing the distributions of the nodes of the graph.
#' @param hypothesis Node to predict.
#' @param path A valid path to save the plots
#' @param directory.name The name of the directory that will be created to save the plots.
#' @return An animated gif that represents the procedure's workflow.
#' @import igraph
#' @examples
#' # load a data set
#' data <- ex1.dag.data
#'
#' # define the distributions of the node
#' mydists <- list(b1="binomial",
#' p1="poisson",
#' g1="gaussian",
#' b2="binomial",
#' p2="poisson",
#' b3="binomial",
#' g2="gaussian",
#' b4="binomial",
#' b5="binomial",
#' g3="gaussian")
#'
#' # infer the graph using ABN
#' max.par <- 4 # set the same max parents for all nodes
#' mycache <- buildScoreCache(data.df = data,
#'                           data.dists = mydists,
#'                           method = "bayes",max.parents = max.par)
#' mp.dag <- mostProbable(score.cache = mycache)
#' dag <- mp.dag$dag
#'
#' hypothesis <- "g2"
#' plot_workflow(dag, mydists, hypothesis)
#' @export

plot_workflow <- function(dag, dists, hypothesis, path = NULL, directory.name =NULL){
  path <- createDirectory(path = path, directory.name = directory.name)

  graph <- graph_from_adjacency_matrix(t(dag))
  node_order <- names(topo_sort(graph, mode="out"))
  node_max <- which(node_order == hypothesis)

  for (i in (1:length(node_order))){
    plots <- plot_Abn(dag, dists, node = node_order[i], "up")
    png(paste0(path,"/graph",i,".png"), width = 800, height = 600)
    renderGraph(plots)
    title(main = paste0("Step ",i,", node ",node_order[i],", up"), cex.main = 1.5, font.main = 2)
    dev.off()
  }

  counter <- length(node_order)
  for (i in (length(node_order):node_max)){
    counter <- counter + 1
    plots <- plot_Abn(dag, dists, node = node_order[i], "down")
    png(paste0(path,"/graph",counter,".png"), width = 800, height = 600)
    renderGraph(plots)
    title(main = paste0("Step ",counter,", node ",node_order[i],", down"), cex.main = 1.5, font.main = 2)
    dev.off()
  }

  createAnimation(path)
}

#' Create a directory
#'
#' This function creates a directory situated in a particular path
#'
#' @param path A valid path.
#' @param directory.name The name of the directory that will be created in this path.
#' @return A path corresponding to this created repository
#' @examples
#' path <- createDirectory(path = NULL, directory.name = "test")
#' @export
createDirectory <- function(path = NULL, directory.name = NULL){
  if (is.null(path)){
    path <- getwd()
  }
  if (is.null(directory.name)){
    if (file.exists("graphs")){
      warning("No directory name was provided, we will use the repository graphs that is already existing. Consider cleaning it before running this code.")
    } else {
      dir.create(file.path(path, "graphs"))
    }
    path <- paste0(path,"/graphs")
  } else {
    if (file.exists(directory.name)){
      warning("The repository already exists. Consider cleaning it before running this code.")
    } else {
      dir.create(file.path(path,directory.name))
    }
    path <- paste0(path,"/",directory.name)
  }
  return(path)
}

#' Create an animated gif
#'
#' This function creates an animated gif file based on png files
#'
#' @param path A valid path to a repository that contains the png to combine in a gif.
#' @return None
#' @importFrom magick image_join
#' @importFrom magick image_read
#' @importFrom magick image_animate
#' @importFrom magick image_write
#' @importFrom purrr map
#' @export
createAnimation <- function(path = NULL){
  png_files <- list.files(path,
                          pattern = "\\.png$",
                          recursive = FALSE,
                          all.files = FALSE,
                          full.names = TRUE)
  png_files <- png_files[order(as.numeric(gsub("\\D", "", png_files)))]

  file.name <- "AnimatedGraph"

  png_files %>%
    map(image_read) %>% # reads each path file
    image_join() %>% # joins image
    image_animate(fps = 1) %>% # animates
    image_write(paste0(path,"/",file.name,".gif"))
}
