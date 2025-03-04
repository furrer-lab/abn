#' Performing inference with ABN
#'
#' Main function to predict the distribution of one node in a fitted ABN graph
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
#' @param dag An adjacency matrix (can be the output of the function mostProbable()).
#' @param myfit Parameters of the network (can be the output of the function fitAbn()).
#' @param hypothesis Node to predict. 
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param plot TRUE/FALSE to indicate if the predicted distribution has to be plotted.
#' @return A list containing the predicted distribution of the hypothesis and predicted distributions of the upstream nodes.
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
#' predictions$predictions_hypothesis # the predicted distribution of g2
#' @export
#' 
predictABN <- function(data, mydists, dag, myfit, hypothesis, evidence, plot = FALSE){
  
  # some checks
  if (ncol(data) != length(mydists) || ncol(data) != ncol(dag) || length(mydists) != ncol(dag)){
    stop("The number of nodes/variables do not correspond.")
  }
  if (!all(colnames(data) %in% names(mydists)) || !all(colnames(data) %in% colnames(dag)) || !all(names(mydists) %in% colnames(dag))){
    stop("The names of the nodes/variables in data, mydists and dag do not correspond.")
  }
  graph <- graph_from_adjacency_matrix(t(dag))
  node_order <- names(topo_sort(graph, mode="out"))

  # Step 0: check the evidence
  evidence <- check_evidence(data, mydists, hypothesis,evidence)

  # Step 1: from top to bottom
  predictions <- list()
  predictions_names <- c()
  for (i in (1:length(node_order))){
    prediction <- predict_node_from_parent(data, mydists, graph, myfit, node = node_order[i], evidence, predictions)
    predictions <- c(predictions, list(prediction))
    predictions_names <- c(predictions_names, node_order[i])
    names(predictions) <- predictions_names
  }

  # Step 2: from bottom to top
  node_max <- which(node_order == hypothesis)
  for (i in (length(node_order):node_max)){
    prediction <- predict_node_from_children(data, mydists, graph, myfit, node = node_order[i], evidence, predictions)
    predictions[[i]] <- prediction
  }
  
  # Step 3: Plot the workflow and the posterior distribution
  if (plot==TRUE){
    path <- createDirectory(path = path, directory.name = directory.name)
    
    for (i in (1:length(node_order))){
      plots <- plot_Abn(dag, mydists, node = node_order[i], "up")
      png(paste0(path,"/graph",i,".png"), width = 800, height = 600)  
      renderGraph(plots)
      title(main = paste0("Step ",i,", node ",node_order[i],", up"), cex.main = 1.5, font.main = 2)
      dev.off()  
    }
    
    counter <- length(node_order)
    for (i in (length(node_order):node_max)){
      counter <- counter + 1
      plots <- plot_Abn(dag, mydists, node = node_order[i], "down")
      png(paste0(path,"/graph",counter,".png"), width = 800, height = 600) 
      renderGraph(plots)
      title(main = paste0("Step ",counter,", node ",node_order[i],", down"), cex.main = 1.5, font.main = 2)
      dev.off()  
    }
    
    createAnimation(path)
    
    g <- plotPosteriorDistrib(predictions,hypothesis,mydists,path)
    print(g)
  }
  return(list(prediction_hypothesis = predictions[[hypothesis]], predictions = predictions))
}

#' Checking the evidence
#'
#' Check the format of the evidence
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
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
check_evidence <- function(data, mydists, hypothesis, evidence){
  if (length(evidence)>0){
    # at least one evidence
    if (hypothesis %in% names(evidence)){
      print(paste0("Node ",hypothesis," is already known (evidence)."))
      evidence <- evidence[hypothesis]
    } else {
      evidence.to.remove <- c()
      for (i in (1:length(evidence))){
        dist.evidence <- mydists[[names(evidence)[i]]]
        if (dist.evidence == "binomial"){
          if (! evidence[[i]] %in% levels(data[[names(evidence[i])]])){
            warning(paste0("Evidence ",names(evidence[i])," does not have an expected value. It should be either ",levels(data[[names(evidence[i])]])[1]," or ",levels(data[[names(evidence[i])]])[2],". It will be discarded."))
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
      if (length(which(mydists[names(evidence)]=="binomial"))>0){
        # at least one evidence is a binomial node
        for (i in (1:length(which(mydists[names(evidence)]=="binomial")))){
          node <- names(which(mydists[names(evidence)]=="binomial"))[i]
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

#' Performing inference with ABN (upstream inference)
#'
#' Main function to predict the distribution of one node given its parents only
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
#' @param graph A dag (igraph object).
#' @param myfit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict. 
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param predictions The estimated predictions of the upstream nodes.
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
#' graph <- graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#' 
#' node <- "g2"
#' 
#' evidence <- list("p1" = 3, "g3" = 5)
#' 
#' predictions <- predict_node_from_parent(data, mydists, graph, myfit, node, evidence)
#' predictions  # the predicted distribution of g2
#' @export
#' 
predict_node_from_parent <- function(data, mydists, graph, myfit, node, evidence, predictions = evidence){
  parents <- find_parents(graph, node)
  if (mydists[[node]]=="poisson"){
    results <- predict_node_from_parent_poisson(data, mydists, myfit, node, evidence, parents, predictions)
  } else if (mydists[[node]]=="gaussian"){
    results <- predict_node_from_parent_gaussian(data, mydists, myfit, node, evidence, parents, predictions)
  } else if (mydists[[node]]=="binomial"){
    results <- predict_node_from_parent_binomial(data, mydists, myfit, node, evidence, parents, predictions)
  }
  return(results)
}

#' Performing inference with ABN (downstream inference)
#'
#' Main function to predict the distribution of one node given its children and its children's parents
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
#' @param graph A dag (igraph object).
#' @param myfit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict. 
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param predictions The estimated predictions of the upstream nodes.
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
#' graph <- graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#' 
#' node <- "g2"
#' 
#' evidence <- list("p1" = 3, "b5" = "y")
#' 
#' predictions <- predict_node_from_children(data, mydists, graph, myfit, node, evidence)
#' predictions  # the predicted distribution of g2
#' @export
#'
predict_node_from_children <- function(data, mydists, graph, myfit, node, evidence, predictions = evidence){
  if (node %in% names(evidence)){
    predictions[[node]]
  } else {
    children <- find_children(graph, node)

    if (length(children)==0){
      # no children
      results <- predictions[[node]]
      return(results)
    } else {
      Results <- list()
      for (j in (1:length(children))){
        child <- children[j]
        parents <- find_parents(graph,child)

        if (mydists[[child]] == "gaussian"){
          results <- predict_node_from_children_gaussian(data, mydists, myfit, node, evidence, child, parents, predictions)
        } else if (mydists[[child]] == "poisson"){
          results <- predict_node_from_children_poisson(data, mydists, myfit, node, evidence, child, parents, predictions)
        } else if (mydists[[child]]=="binomial"){
          results <- predict_node_from_children_binomial(data, mydists, myfit, node, evidence, child, parents, predictions)
        }
        Results <- c(Results,list(results))
      }
      if (mydists[[node]] == "gaussian"){
        results_mean <- mean(sapply(Results,function(l){
          l[1]
        }))
        results_variance <- (1/length(children)) * mean(sapply(Results,function(l){
          l[2]
        }))
        return(c(results_mean,results_variance))
      } else if (mydists[[node]]== "binomial"){
        results_mean <- mean(sapply(Results,function(l){
          l[2]
        }))
        return(c(1-results_mean,results_mean))
      } else {
        results_mean <- mean(sapply(Results,function(l){
          l[1]
        }))
        return(results_mean)
      }
    }
  }
}

#' Performing inference with ABN (upstream inference)
#'
#' Main function to predict the distribution of a Poisson node given its parents only
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
#' @param myfit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict. 
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param parents The parents of the node to predict.
#' @param predictions The estimated predictions of the upstream nodes.
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
#' graph <- graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#' 
#' node <- "p1"
#' parents <- find_parents(graph,node)
#' 
#' evidence <- list("g1" = 4, "g2" = 3, "g3" = 1, b3 = "y")
#' 
#' predictions <- predict_node_from_parent_poisson(data, mydists, myfit, node, evidence, parents)
#' predictions  # the predicted distribution of p1
#' @export
#' 
predict_node_from_parent_poisson <- function(data, mydists, myfit, node, evidence, parents, predictions = evidence){
  if (node %in% names(evidence)){
    # node is an evidence
    node_hat <- evidence[[node]]
  } else {
    if (length(parents)==0){
      # no parents
      node_hat <- predict_root(data, mydists, node)
    } else {
      eq <- myfit[[node]]
      names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])

      bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
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
            proba_tmp <- c(0,0)
            names(proba_tmp) <- levels(data[[bin.nodes.evidence[i]]])
            proba_tmp[grep(probabilities[[bin.nodes.evidence[i]]],names(proba_tmp))] <- 1
            probabilities[[bin.nodes.evidence[i]]] <- proba_tmp
          }
        }

        combinations <- expand.grid(rep(list(c(0, 1)), length(probabilities)))
        names(combinations) <- bin.nodes
        combinations <- as.matrix(combinations)
        if (length(bin.nodes)>1){
          combinations_tmp <- combinations %*% diag(eq[bin.nodes])
        } else {
          combinations_tmp <- combinations %*% eq[bin.nodes]
        }

        proba_cond_values <- apply(combinations_tmp, 1, function(b_vals) {
          exp(continuous_part + sum(b_vals))
        })

        combination_probabilities <- apply(combinations, 1, function(b_vals) {
          prod(sapply(1:length(b_vals), function(i) probabilities[[i]][b_vals[i] + 1]))
        })

        node_hat <- sum(proba_cond_values * combination_probabilities)
      } else {
        node_hat <- exp(continuous_part)
      }
    }
  }
  return(node_hat = node_hat)
}

#' Performing inference with ABN (upstream inference)
#'
#' Main function to predict the distribution of a Gaussian node given its parents only
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
#' @param myfit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict. 
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param parents The parents of the node to predict.
#' @param predictions The estimated predictions of the upstream nodes.
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
#' graph <- graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#' 
#' node <- "g2"
#' parents <- find_parents(graph,node)
#' 
#' evidence <- list("g3" = 1, b3 = "y")
#' 
#' predictions <- predict_node_from_parent_gaussian(data, mydists, myfit, node, evidence, parents)
#' predictions  # the predicted distribution of g2
#' @export
#' 
predict_node_from_parent_gaussian <- function(data, mydists, myfit, node, evidence, parents, predictions = evidence){
  if (node %in% names(evidence)){
    # node is an evidence
    node_hat <- c(evidence[[node]],var(data[[node]]))
  } else {
    if (length(parents)==0){
      # no parents
      node_hat <- predict_root(data, mydists, node)
    } else {
      eq <- myfit[[node]]
      names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])
      bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
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

      # variance estimate
      residuals <- data[node] - (data.matrix(data[,other.nodes]) %*% eq[other.nodes])

      if (length(bin.nodes)>0){
        probabilities <- predictions[bin.nodes]

        bin.nodes.evidence <- intersect(names(evidence),bin.nodes)

        if (length(bin.nodes.evidence)>0){
          # at least one bin nodes is an evidence
          for (i in (1:length(bin.nodes.evidence))){
            proba_tmp <- c(0,0)
            names(proba_tmp) <- levels(data[[bin.nodes.evidence[i]]])
            proba_tmp[grep(probabilities[[bin.nodes.evidence[i]]],names(proba_tmp))] <- 1
            probabilities[[bin.nodes.evidence[i]]] <- proba_tmp
          }
        }

        combinations <- expand.grid(rep(list(c(0, 1)), length(probabilities)))
        names(combinations) <- bin.nodes
        combinations <- as.matrix(combinations)
        if (length(bin.nodes)>1){
          combinations_tmp <- combinations %*% diag(eq[bin.nodes])
        } else {
          combinations_tmp <- combinations %*% eq[bin.nodes]
        }

        proba_cond_values <- apply(combinations_tmp, 1, function(b_vals) {
          continuous_part + sum(b_vals)
        })

        combination_probabilities <- apply(combinations, 1, function(b_vals) {
          prod(sapply(1:length(b_vals), function(i) probabilities[[i]][b_vals[i] + 1]))
        })

        if (length(bin.nodes)>1){
          data_bin <- sapply(data[bin.nodes],function(x){
            ifelse(x == levels(x)[2], 1, 0)
          })
          if (length(bin.nodes.evidence)==0){
            residuals <- residuals - data.matrix(data_bin[,bin.nodes]) %*% eq[bin.nodes]
          } else if (length(setdiff(bin.nodes,bin.nodes.evidence))>1){
            residuals <- residuals - data.matrix(data_bin[,setdiff(bin.nodes,bin.nodes.evidence)]) %*% eq[setdiff(bin.nodes,bin.nodes.evidence)]
            for (i in (1:length(bin.nodes.evidence))){
              residuals <- residuals - as.numeric(c(0,1)[grep(evidence[[bin.nodes.evidence[i]]],levels(data[[bin.nodes.evidence[i]]]))]) * eq[bin.nodes]
            }
          } else if (length(setdiff(bin.nodes,bin.nodes.evidence))==0) {
            for (i in (1:length(bin.nodes.evidence))){
              residuals <- residuals - as.numeric(c(0,1)[grep(evidence[[bin.nodes.evidence[i]]],levels(data[[bin.nodes.evidence[i]]]))]) * eq[bin.nodes]
            }
          } else if (length(setdiff(bin.nodes,bin.nodes.evidence))==0) {
            for (i in (1:length(bin.nodes.evidence))){
              residuals <- residuals - as.numeric(c(0,1)[grep(evidence[[bin.nodes.evidence[i]]],levels(data[[bin.nodes.evidence[i]]]))]) * eq[bin.nodes]
            }
          } else {
            residuals <- residuals - data_bin * eq[setdiff(bin.nodes,bin.nodes.evidence)]
            for (i in (1:length(bin.nodes.evidence))){
              residuals <- residuals - as.numeric(c(0,1)[grep(evidence[[bin.nodes.evidence[i]]],levels(data[[bin.nodes.evidence[i]]]))]) * eq[bin.nodes]
            }
          }
        } else {
          data_bin  <- ifelse(data[[bin.nodes]] == levels(data[[bin.nodes]])[2], 1, 0)
          if (length(bin.nodes.evidence)==0){
            residuals <- residuals - data_bin * eq[bin.nodes]
          } else {
            residuals <- residuals - as.numeric(c(0,1)[grep(evidence[[bin.nodes.evidence]],levels(data[[bin.nodes.evidence]]))]) * eq[bin.nodes]
          }
        }
        node_hat <- c(sum(proba_cond_values * combination_probabilities),sum(residuals^2)/(nrow(residuals)-length(parents)))
      } else {
        node_hat <- c(continuous_part,sum(residuals^2)/(nrow(residuals)-length(parents)))
      }
    }
  }
  return(node_hat = node_hat)
}

#' Performing inference with ABN (upstream inference)
#'
#' Main function to predict the distribution of a binomial node given its parents only
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
#' @param myfit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict. 
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param parents The parents of the node to predict.
#' @param predictions The estimated predictions of the upstream nodes.
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
#' graph <- graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#' 
#' node <- "b2"
#' parents <- find_parents(graph,node)
#' 
#' evidence <- list("g2" = 4, "b1" = "y", "b3" = "y", "b4" = "y")
#' 
#' predictions <- predict_node_from_parent_binomial(data, mydists, myfit, node, evidence, parents)
#' predictions  # the predicted distribution of b2
#' @export
#' 
predict_node_from_parent_binomial <- function(data, mydists, myfit, node, evidence, parents, predictions = evidence){
  if (node %in% names(evidence)){
    # node is an evidence
    node_hat <- evidence[[node]]
  } else {
    if (length(parents)==0){
      # no parents
      node_hat <- predict_root(data, mydists, node)
    } else {
      eq <- unlist(myfit[node])
      names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])
      bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
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
            proba_tmp <- c(0,0)
            names(proba_tmp) <- levels(data[[bin.nodes.evidence[i]]])
            proba_tmp[grep(probabilities[[bin.nodes.evidence[i]]],names(proba_tmp))] <- 1
            probabilities[[bin.nodes.evidence[i]]] <- proba_tmp
          }
        }

        combinations <- expand.grid(rep(list(c(0, 1)), length(probabilities)))
        names(combinations) <- bin.nodes
        combinations <- as.matrix(combinations)
        if (length(bin.nodes)>1){
          combinations_tmp <- combinations %*% diag(eq[bin.nodes])
        } else {
          combinations_tmp <- combinations %*% eq[bin.nodes]
        }

        proba_cond_values <- apply(combinations_tmp, 1, function(b_vals) {
          1 / (1 + exp(continuous_part + sum(b_vals)))
        })
        combination_probabilities <- apply(combinations, 1, function(b_vals) {
          prod(sapply(1:length(b_vals), function(i) probabilities[[i]][b_vals[i] + 1]))
        })

        node_hat <- sum(proba_cond_values * combination_probabilities)
      } else {
        node_hat <- 1/(1+exp(continuous_part))
      }

      node_hat <- c(node_hat,1-node_hat)
      names(node_hat) <- levels(data[[node]])
    }
  }
  return(node_hat)
}

#' Performing inference with ABN (downstream inference)
#'
#' Main function to predict the distribution of a Gaussian node given its children and its children's parents
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
#' @param myfit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict. 
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param children The children of the node to predict.
#' @param parents The children's parents of the node to predict.
#' @param predictions The estimated predictions of the upstream nodes.
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
#' graph <- graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#' 
#' node <- "g2"
#' children <- find_children(graph, node)[1]
#' parents <- find_parents(graph, children)
#' evidence <- list("p1" = 3, "b5" = "y")
#' 
#' predictions <- predict_node_from_children_gaussian(data, mydists, myfit, node, evidence, children, parents)
#' predictions  # the predicted distribution of g2
#' @export
#'
predict_node_from_children_gaussian <- function(data, mydists, myfit, node, evidence, children, parents, predictions){

  if (mydists[[node]] == "binomial"){
    p_prior <- predictions[[node]]

    eq <- myfit[[children]]
    names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])

    bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
    bin.nodes <- setdiff(bin.nodes,node)

    if (length(bin.nodes)>0){
      other.nodes <- parents[-which(parents %in% bin.nodes)]
      other.nodes <- setdiff(other.nodes,node)
    } else {
      other.nodes <- parents
      other.nodes <- setdiff(other.nodes,node)
    }

    predictions_tmp <- predictions[other.nodes]
    predictions_tmp <- lapply(predictions_tmp,function(l){
      l[1]
    })
    continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
    names(continuous_part) <- c()

    if (length(bin.nodes)>0){
      bin.nodes.evidence <- intersect(names(evidence),bin.nodes)

      if (length(bin.nodes.evidence)>0){
        # at least one bin nodes is an evidence
        predictions_tmp <- predictions[bin.nodes.evidence]
        predictions_tmp <- lapply(predictions_tmp,function(l){
          as.numeric(predictions_tmp[[1]])-1
        })
        continuous_part <- continuous_part + sum(eq[bin.nodes.evidence]*unlist(predictions_tmp))

        if (length(bin.nodes.evidence)==length(bin.nodes)){
          # all bin nodes are evidence
          numerator <-  function(x){
            L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part)  * prior_binomial(x, p_prior[2])
          }
          denominator <- numerator(0) + numerator(1)
          results <- c(numerator(0) / denominator,numerator(1) / denominator)
        } else {
          # at least one bin node is not an evidence
          bin.nodes_tmp <- setdiff(bin.nodes,bin.nodes.evidence)
          results_tmp <- c()

          # Generate all combinations of 0 and 1 for binary nodes
          combinations <- expand.grid(rep(list(0:1), length(bin.nodes_tmp)))

          for (i in 1:nrow(combinations)){
            combination_tmp <- as.numeric(combinations[i,])

            # Compute the marginal probability for this combination
            probas <- 1
            for (idx in (1:length(bin.nodes_tmp))) {
              probas <- probas * predictions[[bin.nodes_tmp[idx]]][combination_tmp[idx] + 1]
            }

            # Update the continuous part with contributions from all binary nodes
            continuous_part_tmp <- continuous_part
            for (idx in (1:length(bin.nodes_tmp))) {
              continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes_tmp[idx]] * combination_tmp[idx]
            }

            numerator <-  function(x){
              L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp)  * prior_binomial(x, p_prior[2])
            }
            denominator <- numerator(0) + numerator(1)
            results_tmp <- c(results_tmp,probas*numerator(0) / denominator)
          }
          results <- c(sum(results_tmp),1-sum(results_tmp))
        }
      } else {
        results_tmp <- c()
        # Generate all combinations of 0 and 1 for binary nodes
        combinations <- expand.grid(rep(list(0:1), length(bin.nodes)))

        for (i in 1:nrow(combinations)){
          combination_tmp <- as.numeric(combinations[i,])

          # Compute the marginal probability for this combination
          probas <- 1
          for (idx in (1:length(bin.nodes))) {
            probas <- probas * predictions[[bin.nodes[idx]]][combination_tmp[idx] + 1]
          }

          # Update the continuous part with contributions from all binary nodes
          continuous_part_tmp <- continuous_part
          for (idx in (1:length(bin.nodes))) {
            continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes[idx]] * combination_tmp[idx]
          }

          numerator <-  function(x){
            L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp)  * prior_binomial(x, p_prior[2])
          }
          denominator <- numerator(0) + numerator(1)
          results_tmp <- c(results_tmp,probas*numerator(0) / denominator)
        }
        results <- c(sum(results_tmp),1-sum(results_tmp))
      }
    } else {
      numerator <-  function(x){
        L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part)  * prior_binomial(x, p_prior[2])
      }
      denominator <- numerator(0) + numerator(1)
      results <- c(numerator(0) / denominator,numerator(1) / denominator)
    }
  } else if (mydists[[node]]=="gaussian"){
    mu_prior <- predictions[[node]][1]
    sigma_prior <- predictions[[node]][2]

    eq <- myfit[[children]]
    names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])

    bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
    bin.nodes <- setdiff(bin.nodes,node)

    if (length(bin.nodes)>0){
      other.nodes <- parents[-which(parents %in% bin.nodes)]
      other.nodes <- setdiff(other.nodes,node)
    } else {
      other.nodes <- parents
      other.nodes <- setdiff(other.nodes,node)
    }

    predictions_tmp <- predictions[other.nodes]
    predictions_tmp <- lapply(predictions_tmp,function(l){
      l[1]
    })
    continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
    names(continuous_part) <- c()

    if (length(bin.nodes)>0){
      bin.nodes.evidence <- intersect(names(evidence),bin.nodes)
      if (length(bin.nodes.evidence)>0){
        # at least one bin nodes is an evidence
        predictions_tmp <- predictions[bin.nodes.evidence]
        predictions_tmp <- lapply(predictions_tmp,function(l){
          as.numeric(predictions_tmp[[1]])-1
        })
        continuous_part <- continuous_part + sum(eq[bin.nodes.evidence]*unlist(predictions_tmp))

        if (length(bin.nodes.evidence)==length(bin.nodes)){
          # all bin nodes are evidence
          denominator <- integrate(function(x) L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -Inf, upper = Inf)$value
          numerator <- integrate(function(x) x*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part) * prior_gaussian(x, mu_prior, sigma_prior),
                                 lower = -Inf, upper = Inf)$value
          numerator2 <- integrate(function(x) x^2*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -Inf, upper = Inf)$value
          results <- c(numerator / denominator, numerator2 / denominator - (numerator / denominator)^2)
        } else {
          # at least one bin node is not an evidence
          bin.nodes_tmp <- setdiff(bin.nodes,bin.nodes.evidence)
          results_tmp <- c()
          results_tmp2 <- c()
          
          # Generate all combinations of 0 and 1 for binary nodes
          combinations <- expand.grid(rep(list(0:1), length(bin.nodes_tmp)))

          for (i in 1:nrow(combinations)){
            combination_tmp <- as.numeric(combinations[i,])

            # Compute the marginal probability for this combination
            probas <- 1
            for (idx in (1:length(bin.nodes_tmp))) {
              probas <- probas * predictions[[bin.nodes_tmp[idx]]][combination_tmp[idx] + 1]
            }

            # Update the continuous part with contributions from all binary nodes
            continuous_part_tmp <- continuous_part
            for (idx in (1:length(bin.nodes_tmp))) {
              continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes_tmp[idx]] * combination_tmp[idx]
            }

            denominator <- integrate(function(x) L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior),
                                     lower = -Inf, upper = Inf)$value
            numerator <- integrate(function(x) x*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -Inf, upper = Inf)$value
            numerator2 <- integrate(function(x) x^2*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior),
                                    lower = -Inf, upper = Inf)$value
            results_tmp <- c(results_tmp,probas*(numerator / denominator))
            results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator))
          }
          results <- c(sum(results_tmp),sum(results_tmp2)- (sum(results_tmp))^2)
        }
      } else {
        results_tmp <- c()
        results_tmp2 <- c()
        # Generate all combinations of 0 and 1 for binary nodes
        combinations <- expand.grid(rep(list(0:1), length(bin.nodes)))

        for (i in 1:nrow(combinations)){
          combination_tmp <- as.numeric(combinations[i,])

          # Compute the marginal probability for this combination
          probas <- 1
          for (idx in (1:length(bin.nodes))) {
            probas <- probas * predictions[[bin.nodes[idx]]][combination_tmp[idx] + 1]
          }

          # Update the continuous part with contributions from all binary nodes
          continuous_part_tmp <- continuous_part
          for (idx in (1:length(bin.nodes))) {
            continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes[idx]] * combination_tmp[idx]
          }

          denominator <- integrate(function(x) L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -Inf, upper = Inf)$value
          numerator <- integrate(function(x) x*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior),
                                 lower = -Inf, upper = Inf)$value
          numerator2 <- integrate(function(x) x^2*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior),
                                  lower = -Inf, upper = Inf)$value
          results_tmp <- c(results_tmp,probas*(numerator / denominator))
          results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator))
        }
        results <- c(sum(results_tmp),sum(results_tmp2)- (sum(results_tmp))^2)
      }
    } else {
      denominator <- integrate(function(x) L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part) * prior_gaussian(x, mu_prior, sigma_prior),
                               lower = -Inf, upper = Inf)$value
      numerator <- integrate(function(x) x*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part) * prior_gaussian(x, mu_prior, sigma_prior),
                             lower = -Inf, upper = Inf)$value
      numerator2 <- integrate(function(x) x^2*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part) * prior_gaussian(x, mu_prior, sigma_prior),
                              lower = -Inf, upper = Inf)$value
      results <- c(numerator / denominator, numerator2 / denominator - (numerator / denominator)^2)
    }
  } else {
    lambda_prior <- predictions[[node]]

    eq <- myfit[[children]]
    names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])

    bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
    bin.nodes <- setdiff(bin.nodes,node)

    if (length(bin.nodes)>0){
      other.nodes <- parents[-which(parents %in% bin.nodes)]
      other.nodes <- setdiff(other.nodes,node)
    } else {
      other.nodes <- parents
      other.nodes <- setdiff(other.nodes,node)
    }

    predictions_tmp <- predictions[other.nodes]
    predictions_tmp <- lapply(predictions_tmp,function(l){
      l[1]
    })
    continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
    names(continuous_part) <- c()

    if (length(bin.nodes)>0){
      bin.nodes.evidence <- intersect(names(evidence),bin.nodes)
      if (length(bin.nodes.evidence)>0){
        # at least one bin nodes is an evidence
        predictions_tmp <- predictions[bin.nodes.evidence]
        predictions_tmp <- lapply(predictions_tmp,function(l){
          as.numeric(predictions_tmp[[1]])-1
        })
        continuous_part <- continuous_part + sum(eq[bin.nodes.evidence]*unlist(predictions_tmp))

        if (length(bin.nodes.evidence)==length(bin.nodes)){
          # all bin nodes are evidence
          max_x <- max(1000,4*lambda_prior)
          denominator <- sum(sapply(0:max_x,function(x) L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part) * prior_poisson(x, lambda_prior)))
          numerator <- sum(sapply(0:max_x,function(x) x*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part) * prior_poisson(x, lambda_prior)))
          results <- numerator / denominator
        } else {
          # at least one bin node is not an evidence
          bin.nodes_tmp <- setdiff(bin.nodes,bin.nodes.evidence)
          results_tmp <- c()

          # Generate all combinations of 0 and 1 for binary nodes
          combinations <- expand.grid(rep(list(0:1), length(bin.nodes_tmp)))

          for (i in 1:nrow(combinations)){
            combination_tmp <- as.numeric(combinations[i,])

            # Compute the marginal probability for this combination
            probas <- 1
            for (idx in (1:length(bin.nodes_tmp))) {
              probas <- probas * predictions[[bin.nodes_tmp[idx]]][combination_tmp[idx] + 1]
            }

            # Update the continuous part with contributions from all binary nodes
            continuous_part_tmp <- continuous_part
            for (idx in (1:length(bin.nodes_tmp))) {
              continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes_tmp[idx]] * combination_tmp[idx]
            }
            max_x <- max(1000,4*lambda_prior)
            denominator <- sum(sapply(0:max_x,function(x) L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp) * prior_poisson(x, lambda_prior)))
            numerator <- sum(sapply(0:max_x,function(x) x*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp) * prior_poisson(x, lambda_prior)))
            results_tmp <- c(results_tmp,probas*(numerator / denominator))
          }
          results <- sum(results_tmp)
        }
      } else {
        results_tmp <- c()
        # Generate all combinations of 0 and 1 for binary nodes
        combinations <- expand.grid(rep(list(0:1), length(bin.nodes)))

        for (i in 1:nrow(combinations)){
          combination_tmp <- as.numeric(combinations[i,])

          # Compute the marginal probability for this combination
          probas <- 1
          for (idx in (1:length(bin.nodes))) {
            probas <- probas * predictions[[bin.nodes[idx]]][combination_tmp[idx] + 1]
          }

          # Update the continuous part with contributions from all binary nodes
          continuous_part_tmp <- continuous_part
          for (idx in (1:length(bin.nodes))) {
            continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes[idx]] * combination_tmp[idx]
          }

          max_x <- max(1000,4*lambda_prior)
          denominator <- sum(sapply(0:max_x,function(x) L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp) * prior_poisson(x, lambda_prior)))
          numerator <- sum(sapply(0:max_x, function(x) x*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part_tmp) * prior_poisson(x, lambda_prior)))
          results_tmp <- c(results_tmp,probas*(numerator / denominator))
        }
        results <- sum(results_tmp)
      }
    } else {
      max_x <- max(1000,4*lambda_prior)
      denominator <- sum(sapply(0:max_x,function(x) L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part) * prior_poisson(x, lambda_prior)))
      numerator <- sum(sapply(0:max_x,function(x) x*L_gaussian(y = predictions[[children]][[1]], x, coef = eq[[node]], var = predictions[[children]][2],continuous_part) * prior_poisson(x, lambda_prior)))
      results <- numerator / denominator
    }
  }
}

#' Performing inference with ABN (downstream inference)
#'
#' Main function to predict the distribution of a Poisson node given its children and its children's parents
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
#' @param myfit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict. 
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param children The children of the node to predict.
#' @param parents The children's parents of the node to predict.
#' @param predictions The estimated predictions of the upstream nodes.
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
#' graph <- graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#' 
#' node <- "g2"
#' children <- find_children(graph, node)[1]
#' parents <- find_parents(graph, children)
#' evidence <- list("p1" = 3, "b5" = "y")
#' 
#' predictions <- predict_node_from_children_gaussian(data, mydists, myfit, node, evidence, children, parents)
#' predictions  # the predicted distribution of g2
#' @export
#'
predict_node_from_children_poisson <- function(data, mydists, myfit, node, evidence, children, parents, predictions){

  if (mydists[[node]] == "gaussian"){
    mu_prior <- predictions[[node]][1]
    sigma_prior <- predictions[[node]][2]

    eq <- myfit[[children]]
    
    names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])

    bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
    bin.nodes <- setdiff(bin.nodes,node)

    if (length(bin.nodes)>0){
      other.nodes <- parents[-which(parents %in% bin.nodes)]
      other.nodes <- setdiff(other.nodes,node)
    } else {
      other.nodes <- parents
      other.nodes <- setdiff(other.nodes,node)
    }

    predictions_tmp <- predictions[other.nodes]
    predictions_tmp <- lapply(predictions_tmp,function(l){
      l[1]
    })
    continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
    names(continuous_part) <- c()

    if (length(bin.nodes)>0){
      bin.nodes.evidence <- intersect(names(evidence),bin.nodes)
      if (length(bin.nodes.evidence)>0){
        # at least one bin nodes is an evidence
        predictions_tmp <- predictions[bin.nodes.evidence]
        predictions_tmp <- lapply(predictions_tmp,function(l){
          as.numeric(predictions_tmp[[1]])-1
        })
        continuous_part <- continuous_part + sum(eq[bin.nodes.evidence]*unlist(predictions_tmp))

        if (length(bin.nodes.evidence)==length(bin.nodes)){
          # all bin nodes are evidence
          try.denominator <- try(integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -Inf, upper = Inf),TRUE)
          if (length(try.denominator)==1){
            denominator <- integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                     lower = -10, upper = 10)$value
          } else if (try.denominator$value == 0){
            denominator <- integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                     lower = -10, upper = 10)$value
          } else {
            denominator <- try.denominator$value
          }

          try.numerator <- try(integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                 lower = -Inf, upper = Inf),TRUE)
          if (length(try.numerator)==1){
            numerator <- integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -10, upper = 10)$value
          } else if (try.numerator$value == 0){
            numerator <- integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -10, upper = 10)$value
          } else {
            numerator <- try.numerator$value
          }
          try.numerator2 <- try(integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                         lower = -Inf, upper = Inf),TRUE)
          if (length(try.numerator2)==1){
            numerator2 <- integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -10, upper = 10)$value
          } else if (try.numerator2$value == 0){ 
            numerator2 <- integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                    lower = -10, upper = 10)$value
          } else {
            numerator2 <- try.numerator2$value
          }
          results <- c(numerator / denominator, numerator2 / denominator)
        } else {
          # at least one bin node is not an evidence
          bin.nodes_tmp <- setdiff(bin.nodes,bin.nodes.evidence)
          results_tmp <- c()
          results_tmp2 <- c()
          
          # Generate all combinations of 0 and 1 for binary nodes
          combinations <- expand.grid(rep(list(0:1), length(bin.nodes_tmp)))

          for (i in 1:nrow(combinations)){
            combination_tmp <- as.numeric(combinations[i,])

            # Compute the marginal probability for this combination
            probas <- 1
            for (idx in (1:length(bin.nodes_tmp))) {
              probas <- probas * predictions[[bin.nodes_tmp[idx]]][combination_tmp[idx] + 1]
            }

            # Update the continuous part with contributions from all binary nodes
            continuous_part_tmp <- continuous_part
            for (idx in (1:length(bin.nodes_tmp))) {
              continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes_tmp[idx]] * combination_tmp[idx]
            }

            try.denominator <- try(integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                             lower = -Inf, upper = Inf),TRUE)
            if (length(try.denominator)==1){
              denominator <- integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                       lower = -10, upper = 10)$value
            } else if (try.denominator$value == 0){ 
              denominator <- integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                       lower = -10, upper = 10)$value
            } else {
              denominator <- try.denominator$value
            }

            try.numerator <- try(integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                           lower = -Inf, upper = Inf),TRUE)
            if (length(try.numerator)==1){
              numerator <- integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                     lower = -10, upper = 10)$value
            } else if (try.numerator$value == 0){ 
              numerator <- integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                     lower = -10, upper = 10)$value
            } else {
              numerator <- try.numerator$value
            }

            try.numerator2 <- try(integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                           lower = -Inf, upper = Inf),TRUE)
            if (length(try.numerator2)==1){
              numerator2 <- integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                     lower = -10, upper = 10)$value
            } else if (try.numerator2$value==0){
              numerator2 <- integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                      lower = -10, upper = 10)$value
            } else {
              numerator2 <- try.numerator2$value
            }
            results_tmp <- c(results_tmp,probas*(numerator / denominator))
            results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator))
          }
          results <- c(sum(results_tmp),sum(results_tmp2)- (sum(results_tmp))^2)
        }
      } else {
        results_tmp <- c()
        results_tmp2 <- c()
        # Generate all combinations of 0 and 1 for binary nodes
        combinations <- expand.grid(rep(list(0:1), length(bin.nodes)))

        for (i in 1:nrow(combinations)){
          combination_tmp <- as.numeric(combinations[i,])

          # Compute the marginal probability for this combination
          probas <- 1
          for (idx in (1:length(bin.nodes))) {
            probas <- probas * predictions[[bin.nodes[idx]]][combination_tmp[idx] + 1]
          }

          # Update the continuous part with contributions from all binary nodes
          continuous_part_tmp <- continuous_part
          for (idx in (1:length(bin.nodes))) {
            continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes[idx]] * combination_tmp[idx]
          }

          try.denominator <- try(integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                           lower = -Inf, upper = Inf),TRUE)
          if (length(try.denominator)==1){
            denominator <- integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                     lower = -10, upper = 10)$value
          } else if (try.denominator$value == 0){
            denominator <- integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                     lower = -10, upper = 10)$value
          } else {
            denominator <- try.denominator$value
          }

          try.numerator <- try(integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                         lower = -Inf, upper = Inf),TRUE)
          if (length(try.numerator)==1){
            numerator <- integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -10, upper = 10)$value
          } else if (try.numerator$value==0){
            numerator <- integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -10, upper = 10)$value
          } else {
            numerator <- try.numerator$value
          }

          try.numerator2 <- try(integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                         lower = -Inf, upper = Inf),TRUE)
          if (length(try.numerator2)==1){
            numerator2 <- integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -10, upper = 10)$value
          } else if (try.numerator2$value==0){
            try.numerator2 <- try(integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -10, upper = 10), TRUE)
            if (length(try.numerator2)==1){
              numerator2 <- 0
            } else {
              numerator2 <- integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                          lower = -10, upper = 10)$value
            }
          } else {
            numerator2 <- try.numerator2$value
          }
          
          results_tmp <- c(results_tmp,probas*(numerator / denominator))
          results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator))
        }
        results <- c(sum(results_tmp),sum(results_tmp2)- (sum(results_tmp))^2)
        if (length(which(is.na(results)))>0){
          warning(paste0("Numerical issues with the computation of the variance of node ",node,". The expectation and the variance will not be updated."))
          results[1] <- predictions[[node]][1]
          results[2] <- predictions[[node]][2]
        } else if  (results[2]<=1e-2){
          warning(paste0("Numerical issues with the computation of the variance of node ",node,". The expectation and the variance will not be updated."))
          results[1] <- predictions[[node]][1]
          results[2] <- predictions[[node]][2]
        }
        return(results)
      }
    } else {
      try.denominator <- try(integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                       lower = -Inf, upper = Inf),TRUE)
      if (length(try.denominator)==1){
        denominator <- integrate(function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                 lower = -10, upper = 10)$value
      } else {
        denominator <- try.denominator$value
      }

      try.numerator <- try(integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                     lower = -Inf, upper = Inf),TRUE)
      if (length(try.numerator)==1){
        numerator <- integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                               lower = -10, upper = 10)$value
      } else {
        numerator <- try.numerator$value
      }

      try.numerator2 <- try(integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                                     lower = -Inf, upper = Inf),TRUE)
      if (length(try.numerator2)==1){
        numerator2 <- integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_gaussian(x, mu_prior, sigma_prior),
                               lower = -10, upper = 10)$value
      } else {
        numerator2 <- try.numerator2$value
      }
      results <- c(numerator / denominator, numerator2 / denominator - (numerator / denominator)^2)
    }
  } else if (mydists[[node]] == "binomial"){
    p_prior <- predictions[[node]]

    eq <- myfit[[children]]
    names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])

    bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
    bin.nodes <- setdiff(bin.nodes,node)

    if (length(bin.nodes)>0){
      other.nodes <- parents[-which(parents %in% bin.nodes)]
      other.nodes <- setdiff(other.nodes,node)
    } else {
      other.nodes <- parents
      other.nodes <- setdiff(other.nodes,node)
    }

    predictions_tmp <- predictions[other.nodes]
    predictions_tmp <- lapply(predictions_tmp,function(l){
      l[1]
    })
    continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
    names(continuous_part) <- c()

    if (length(bin.nodes)>0){
      bin.nodes.evidence <- intersect(names(evidence),bin.nodes)

      if (length(bin.nodes.evidence)>0){
        # at least one bin nodes is an evidence
        predictions_tmp <- predictions[bin.nodes.evidence]
        predictions_tmp <- lapply(predictions_tmp,function(l){
          as.numeric(predictions_tmp[[1]])-1
        })
        continuous_part <- continuous_part + sum(eq[bin.nodes.evidence]*unlist(predictions_tmp))

        if (length(bin.nodes.evidence)==length(bin.nodes)){
          # all bin nodes are evidence
          numerator <-  function(x){
            exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]],continuous_part))  * prior_binomial(x, p_prior[2])
          }

          denominator <- numerator(0) + numerator(1)
          results <- c(numerator(0) / denominator,numerator(1) / denominator)
        } else {
          # at least one bin node is not an evidence
          bin.nodes_tmp <- setdiff(bin.nodes,bin.nodes.evidence)
          results_tmp <- c()

          # Generate all combinations of 0 and 1 for binary nodes
          combinations <- expand.grid(rep(list(0:1), length(bin.nodes_tmp)))

          for (i in 1:nrow(combinations)){
            combination_tmp <- as.numeric(combinations[i,])

            # Compute the marginal probability for this combination
            probas <- 1
            for (idx in (1:length(bin.nodes_tmp))) {
              probas <- probas * predictions[[bin.nodes_tmp[idx]]][combination_tmp[idx] + 1]
            }

            # Update the continuous part with contributions from all binary nodes
            continuous_part_tmp <- continuous_part
            for (idx in (1:length(bin.nodes_tmp))) {
              continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes_tmp[idx]] * combination_tmp[idx]
            }

            numerator <-  function(x){
              exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]],continuous_part_tmp))  * prior_binomial(x, p_prior[2])
            }
            denominator <- numerator(0) + numerator(1)
            results_tmp <- c(results_tmp,probas*numerator(0) / denominator)
          }
          results <- c(sum(results_tmp),1-sum(results_tmp))
        }
      } else {
        results_tmp <- c()
        # Generate all combinations of 0 and 1 for binary nodes
        combinations <- expand.grid(rep(list(0:1), length(bin.nodes)))

        for (i in 1:nrow(combinations)){
          combination_tmp <- as.numeric(combinations[i,])

          # Compute the marginal probability for this combination
          probas <- 1
          for (idx in (1:length(bin.nodes))) {
            probas <- probas * predictions[[bin.nodes[idx]]][combination_tmp[idx] + 1]
          }

          # Update the continuous part with contributions from all binary nodes
          continuous_part_tmp <- continuous_part
          for (idx in (1:length(bin.nodes))) {
            continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes[idx]] * combination_tmp[idx]
          }

          numerator <-  function(x){
            exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]],continuous_part_tmp))  * prior_binomial(x, p_prior[2])
          }
          denominator <- numerator(0) + numerator(1)
          results_tmp <- c(results_tmp,probas*numerator(0) / denominator)
        }
        results <- c(sum(results_tmp),1-sum(results_tmp))
      }
    } else {
      numerator <-  function(x){
        exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]],continuous_part)) * prior_binomial(x, p_prior[2])
      }
      denominator <- numerator(0) + numerator(1)
      results <- c(numerator(0) / denominator,numerator(1) / denominator)
      if (length(which(is.na(results)))>0){
        warning(paste0("Numerical issues with the computation of the variance of node ",node,". The distribution will not be updated."))
        results <- predictions[[node]]
      }
      return(results)
    }
  } else {
    lambda_prior <- predictions[[node]]

    eq <- myfit[[children]]
    names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])

    bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
    bin.nodes <- setdiff(bin.nodes,node)

    if (length(bin.nodes)>0){
      other.nodes <- parents[-which(parents %in% bin.nodes)]
      other.nodes <- setdiff(other.nodes,node)
    } else {
      other.nodes <- parents
      other.nodes <- setdiff(other.nodes,node)
    }

    predictions_tmp <- predictions[other.nodes]
    predictions_tmp <- lapply(predictions_tmp,function(l){
      l[1]
    })
    continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
    names(continuous_part) <- c()

    if (length(bin.nodes)>0){
      bin.nodes.evidence <- intersect(names(evidence),bin.nodes)
      if (length(bin.nodes.evidence)>0){
        # at least one bin nodes is an evidence
        predictions_tmp <- predictions[bin.nodes.evidence]
        predictions_tmp <- lapply(predictions_tmp,function(l){
          as.numeric(predictions_tmp[[1]])-1
        })
        continuous_part <- continuous_part + sum(eq[bin.nodes.evidence]*unlist(predictions_tmp))

        if (length(bin.nodes.evidence)==length(bin.nodes)){
          # all bin nodes are evidence
          max_x <- max(1000,4*lambda_prior)
          #try.denominator <- try(sum(sapply(0:max_x,function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_poisson(x,lambda_prior))))
          #if (length(try.denominator)==1){
            denominator <- sum(sapply(0:max_x, function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_poisson(x,lambda_prior)))
          #} else if (try.denominator$value == 0){
          #  denominator <- sum(sapply(0:max_x,function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_poisson(x,lambda_prior)))
          #} else {
          #  denominator <- try.denominator$value
          #}

          #try.numerator <- try(sum(sapply(0:max_x,function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_poisson(x,lambda_prior))))
          #if (length(try.numerator)==1){
            numerator <- sum(sapply(0:max_x,function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_poisson(x,lambda_prior)))
          #} else {
          #  numerator <- try.numerator$value
          #}
          results <- numerator / denominator
        } else {
          # at least one bin node is not an evidence
          bin.nodes_tmp <- setdiff(bin.nodes,bin.nodes.evidence)
          results_tmp <- c()

          # Generate all combinations of 0 and 1 for binary nodes
          combinations <- expand.grid(rep(list(0:1), length(bin.nodes_tmp)))

          for (i in 1:nrow(combinations)){
            combination_tmp <- as.numeric(combinations[i,])

            # Compute the marginal probability for this combination
            probas <- 1
            for (idx in (1:length(bin.nodes_tmp))) {
              probas <- probas * predictions[[bin.nodes_tmp[idx]]][combination_tmp[idx] + 1]
            }

            # Update the continuous part with contributions from all binary nodes
            continuous_part_tmp <- continuous_part
            for (idx in (1:length(bin.nodes_tmp))) {
              continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes_tmp[idx]] * combination_tmp[idx]
            }

            max_x <- max(1000,4*lambda_prior)
            #try.denominator <- try(sum(sapply(0:max_x, function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_poisson(x, lambda_prior))),TRUE)
            #if (length(try.denominator)==1){
              denominator <- sum(sapply(0:max_x,function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_poisson(x, lambda_prior)))
            #} else {
            #  denominator <- try.denominator$value
            #}

            #try.numerator <- try(sum(sapply(0:max_x,function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_poisson(x, lambda_prior))),TRUE)
            #if (length(try.numerator)==1){
              numerator <- sum(sapply(0:max_x,function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_poisson(x, lambda_prior)))
            #} else {
            #  numerator <- try.numerator$value
            #}

            results_tmp <- c(results_tmp,probas*(numerator / denominator))
          }
          results <- sum(results_tmp)
        }
      } else {
        results_tmp <- c()
        # Generate all combinations of 0 and 1 for binary nodes
        combinations <- expand.grid(rep(list(0:1), length(bin.nodes)))

        for (i in 1:nrow(combinations)){
          combination_tmp <- as.numeric(combinations[i,])

          # Compute the marginal probability for this combination
          probas <- 1
          for (idx in (1:length(bin.nodes))) {
            probas <- probas * predictions[[bin.nodes[idx]]][combination_tmp[idx] + 1]
          }

          # Update the continuous part with contributions from all binary nodes
          continuous_part_tmp <- continuous_part
          for (idx in (1:length(bin.nodes))) {
            continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes[idx]] * combination_tmp[idx]
          }

          max_x <- max(1000,4*lambda_prior)
          #try.denominator <- try(sum(sapply(0:max_x,function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_poisson(x, lambda_prior))),TRUE)
          #if (length(try.denominator)==1){
            denominator <- sum(sapply(0:max_x,function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_poisson(x, lambda_prior)))
          #} else if (try.denominator$value == 0){
          #  denominator <- sum(sapply(0:max_x,function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_poisson(x, lambda_prior)))
          #} else {
          #  denominator <- try.denominator$value
          #}

          #try.numerator <- try(sum(sapply(0:max_x,function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_poisson(x, lambda_prior))),TRUE)
          #if (length(try.numerator)==1){
            numerator <- sum(sapply(0:max_x,function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_poisson(x, lambda_prior)))
          #} else if (try.numerator$value==0){
          #  numerator <- sum(sapply(0:max_x,function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_poisson(x, lambda_prior)))
          #} else {
          #  numerator <- try.numerator$value
          #}

          results_tmp <- c(results_tmp,probas*(numerator / denominator))
        }
        results <- sum(results_tmp)
      }
    } else {
      max_x <- max(1000,4*lambda_prior)
      #try.denominator <- try(sum(sapply(0:max_x,function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_poisson(x, lambda_prior))),TRUE)
      #if (length(try.denominator)==1){
        denominator <- sum(sapply(0:max_x,function(x) exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_poisson(x, lambda_prior)))
      #} else {
      #  denominator <- try.denominator$value
      #}

      #try.numerator <- try(sum(sapply(0:max_x,function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_poisson(x, lambda_prior))),TRUE)
      #if (length(try.numerator)==1){
        numerator <- sum(sapply(0:max_x,function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part)) * prior_poisson(x, lambda_prior)))
      #} else {
      #  numerator <- try.numerator$value
      #}

      results <- numerator / denominator
    }
  }
}

#' Performing inference with ABN (downstream inference)
#'
#' Main function to predict the distribution of a binomial node given its children and its children's parents
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
#' @param myfit Parameters of the network (can be the output of the function fitAbn()).
#' @param node Temporary node to predict. 
#' @param evidence Known nodes that are used to predict the hypothesis.
#' @param children The children of the node to predict.
#' @param parents The children's parents of the node to predict.
#' @param predictions The estimated predictions of the upstream nodes.
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
#' graph <- graph_from_adjacency_matrix(t(dag))
#'
#' # infer the parameters of the network
#' myfit <- fitAbn(object = mp.dag)
#' myfit <- myfit$modes
#' 
#' node <- "g2"
#' children <- find_children(graph, node)[1]
#' parents <- find_parents(graph, children)
#' evidence <- list("p1" = 3, "b5" = "y")
#' 
#' predictions <- predict_node_from_children_gaussian(data, mydists, myfit, node, evidence, children, parents)
#' predictions  # the predicted distribution of g2
#' @export
#'
predict_node_from_children_binomial <- function(data, mydists, myfit, node, evidence, children, parents, predictions){

  if (mydists[[node]] == "binomial"){
    p_prior <- predictions[[node]]

    eq <- myfit[[children]]
    names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])

    bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
    bin.nodes <- setdiff(bin.nodes,node)

    if (length(bin.nodes)>0){
      other.nodes <- parents[-which(parents %in% bin.nodes)]
      other.nodes <- setdiff(other.nodes,node)
    } else {
      other.nodes <- parents
      other.nodes <- setdiff(other.nodes,node)
    }

    predictions_tmp <- predictions[other.nodes]
    predictions_tmp <- lapply(predictions_tmp,function(l){
      l[1]
    })
    continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
    names(continuous_part) <- c()

    if (children %in% names(evidence)){
      p_pred <- c(0,0)
      names(p_pred) <- levels(data[[children]])
      p_pred[grep(predictions[[children]],names(p_pred))] <- 1
    } else {
      p_pred <- predictions[[children]]
    }
    
    if (length(bin.nodes)>0){
      bin.nodes.evidence <- intersect(names(evidence),bin.nodes)

      if (length(bin.nodes.evidence)>0){
        # at least one bin nodes is an evidence
        predictions_tmp <- predictions[bin.nodes.evidence]
        predictions_tmp <- lapply(predictions_tmp,function(l){
          as.numeric(l[[1]])-1
        })
        continuous_part <- continuous_part + sum(eq[bin.nodes.evidence]*unlist(predictions_tmp))

        if (length(bin.nodes.evidence)==length(bin.nodes)){
          # all bin nodes are evidence
          numerator <-  function(x,y){
            L_binomial(y, x, coef = eq[[node]], continuous_part)  * prior_binomial(x, p_prior[2])
          }
          
          numerator2 <- function(x){
            numerator(x,0)*p_pred[1] + numerator(x,1)*p_pred[2]
          }
          
          denominator <- numerator2(0) + numerator2(1)
          results <- c(numerator2(0) / denominator,numerator2(1) / denominator)
        } else {
          # at least one bin node is not an evidence
          bin.nodes_tmp <- setdiff(bin.nodes,bin.nodes.evidence)
          results_tmp <- c()

          # Generate all combinations of 0 and 1 for binary nodes
          combinations <- expand.grid(rep(list(0:1), length(bin.nodes_tmp)))

          for (i in 1:nrow(combinations)){
            combination_tmp <- as.numeric(combinations[i,])

            # Compute the marginal probability for this combination
            probas <- 1
            for (idx in (1:length(bin.nodes_tmp))) {
              probas <- probas * predictions[[bin.nodes_tmp[idx]]][combination_tmp[idx] + 1]
            }

            # Update the continuous part with contributions from all binary nodes
            continuous_part_tmp <- continuous_part
            for (idx in (1:length(bin.nodes_tmp))) {
              continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes_tmp[idx]] * combination_tmp[idx]
            }

            numerator <-  function(x,y){
              L_binomial(y, x, coef = eq[[node]], continuous_part_tmp)  * prior_binomial(x, p_prior[2])
            }
            
            numerator2 <- function(x){
              numerator(x,0)*p_pred[1] + numerator(x,1)*p_pred[2]
            }

            denominator <- numerator2(0) + numerator2(1)
            results_tmp <- c(results_tmp,probas*numerator2(0) / denominator)
          }
          results <- c(sum(results_tmp),1-sum(results_tmp))
        }
      } else {
        results_tmp <- c()
        # Generate all combinations of 0 and 1 for binary nodes
        combinations <- expand.grid(rep(list(0:1), length(bin.nodes)))

        for (i in 1:nrow(combinations)){
          combination_tmp <- as.numeric(combinations[i,])

          # Compute the marginal probability for this combination
          probas <- 1
          for (idx in (1:length(bin.nodes))) {
            probas <- probas * predictions[[bin.nodes[idx]]][combination_tmp[idx] + 1]
          }

          # Update the continuous part with contributions from all binary nodes
          continuous_part_tmp <- continuous_part
          for (idx in (1:length(bin.nodes))) {
            continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes[idx]] * combination_tmp[idx]
          }

          numerator <-  function(x,y){
            L_binomial(y, x, coef = eq[[node]], continuous_part_tmp)  * prior_binomial(x, p_prior[2])
          }

          numerator2 <- function(x){
            numerator(x,0)*p_pred[1] + numerator(x,1)*p_pred[2]
          }
          
          denominator <- numerator2(0) + numerator2(1)
          results_tmp <- c(results_tmp,probas*numerator2(0) / denominator)
        }
        results <- c(sum(results_tmp),1-sum(results_tmp))
      }
    } else {
      numerator <-  function(x,y){
        L_binomial(y, x, coef = eq[[node]], continuous_part)  * prior_binomial(x, p_prior[2])
      }
      
      numerator2 <- function(x){
        numerator(x,0)*p_pred[1] + numerator(x,1)*p_pred[2]
      }
      
      denominator <- numerator2(0) + numerator2(1)
      results <- c(numerator2(0) / denominator,numerator2(1) / denominator)
    }
  } else if (mydists[[node]]=="gaussian"){
    mu_prior <- predictions[[node]][1]
    sigma_prior <- predictions[[node]][2]

    eq <- myfit[[children]]
    names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])

    bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
    bin.nodes <- setdiff(bin.nodes,node)
    if (length(bin.nodes)>0){
      other.nodes <- parents[-which(parents %in% bin.nodes)]
      other.nodes <- setdiff(other.nodes,node)
    } else {
      other.nodes <- parents
      other.nodes <- setdiff(other.nodes,node)
    }

    predictions_tmp <- predictions[other.nodes]
    predictions_tmp <- lapply(predictions_tmp,function(l){
      l[1]
    })
    continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
    names(continuous_part) <- c()

    if (children %in% names(evidence)){
      p_prior <- c(0,0)
      names(p_prior) <- levels(data[[children]])
      p_prior[grep(predictions[[children]],names(p_prior))] <- 1
    } else {
      p_prior <- predictions[[children]]
    }

    if (length(bin.nodes)>0){
      bin.nodes.evidence <- intersect(names(evidence),bin.nodes)
      if (length(bin.nodes.evidence)>0){
        # at least one bin nodes is an evidence
        predictions_tmp <- predictions[bin.nodes.evidence]
        predictions_tmp <- lapply(predictions_tmp,function(l){
          as.numeric(predictions_tmp[[1]])-1
        })
        continuous_part <- continuous_part + sum(eq[bin.nodes.evidence]*unlist(predictions_tmp))

        if (length(bin.nodes.evidence)==length(bin.nodes)){
          # all bin nodes are evidence
          denominator <- sum(c(integrate(function(x) L_binomial(y = 0, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                                         lower = -Inf, upper = Inf)$value,
                               integrate(function(x) L_binomial(y = 1, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                         lower = -Inf, upper = Inf)$value))
          numerator <- sum(c(integrate(function(x) x*L_binomial(y = 0, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                                       lower = -Inf, upper = Inf)$value,
                             integrate(function(x) x*L_binomial(y = 1, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                       lower = -Inf, upper = Inf)$value))
          numerator2 <- sum(c(integrate(function(x) x^2*L_binomial(y = 0, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                                        lower = -Inf, upper = Inf)$value,
                              integrate(function(x) x^2*L_binomial(y = 1, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                        lower = -Inf, upper = Inf)$value))
          results <- c(numerator / denominator, numerator2 / denominator - (numerator / denominator)^2)
        } else {
          # at least one bin node is not an evidence
          bin.nodes_tmp <- setdiff(bin.nodes,bin.nodes.evidence)
          results_tmp <- c()
          results_tmp2 <- c()
          
          # Generate all combinations of 0 and 1 for binary nodes
          combinations <- expand.grid(rep(list(0:1), length(bin.nodes_tmp)))

          for (i in 1:nrow(combinations)){
            combination_tmp <- as.numeric(combinations[i,])

            # Compute the marginal probability for this combination
            probas <- 1
            for (idx in (1:length(bin.nodes_tmp))) {
              probas <- probas * predictions[[bin.nodes_tmp[idx]]][combination_tmp[idx] + 1]
            }

            # Update the continuous part with contributions from all binary nodes
            continuous_part_tmp <- continuous_part
            for (idx in (1:length(bin.nodes_tmp))) {
              continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes_tmp[idx]] * combination_tmp[idx]
            }

            denominator <- sum(c(integrate(function(x) L_binomial(y = 0, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                                           lower = -Inf, upper = Inf)$value,
                                 integrate(function(x) L_binomial(y = 1, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                           lower = -Inf, upper = Inf)$value))
            numerator <- sum(c(integrate(function(x) x*L_binomial(y = 0, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                                         lower = -Inf, upper = Inf)$value,
                               integrate(function(x) x*L_binomial(y = 1, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                         lower = -Inf, upper = Inf)$value))
            numerator2 <- sum(c(integrate(function(x) x^2*L_binomial(y = 0, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                                          lower = -Inf, upper = Inf)$value,
                                integrate(function(x) x^2*L_binomial(y = 1, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                          lower = -Inf, upper = Inf)$value))
            
            results_tmp <- c(results_tmp,probas*(numerator / denominator))
            results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator))
          }
          results <- c(sum(results_tmp),sum(results_tmp2)- (sum(results_tmp))^2)
        }
      } else {
        results_tmp <- c()
        results_tmp2 <- c()
        # Generate all combinations of 0 and 1 for binary nodes
        combinations <- expand.grid(rep(list(0:1), length(bin.nodes)))

        for (i in 1:nrow(combinations)){
          combination_tmp <- as.numeric(combinations[i,])

          # Compute the marginal probability for this combination
          probas <- 1
          for (idx in (1:length(bin.nodes))) {
            probas <- probas * predictions[[bin.nodes[idx]]][combination_tmp[idx] + 1]
          }

          # Update the continuous part with contributions from all binary nodes
          continuous_part_tmp <- continuous_part
          for (idx in (1:length(bin.nodes))) {
            continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes[idx]] * combination_tmp[idx]
          }

          denominator <- sum(c(integrate(function(x) L_binomial(y = 0, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                                         lower = -Inf, upper = Inf)$value,
                               integrate(function(x) L_binomial(y = 1, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                         lower = -Inf, upper = Inf)$value))
          numerator <- sum(c(integrate(function(x) x*L_binomial(y = 0, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                                       lower = -Inf, upper = Inf)$value,
                             integrate(function(x) x*L_binomial(y = 1, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                       lower = -Inf, upper = Inf)$value))
          numerator2 <- sum(c(integrate(function(x) x^2*L_binomial(y = 0, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                                        lower = -Inf, upper = Inf)$value,
                              integrate(function(x) x^2*L_binomial(y = 1, x, coef = eq[[node]], continuous_part_tmp) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                        lower = -Inf, upper = Inf)$value))
          
          results_tmp <- c(results_tmp,probas*(numerator / denominator))
          results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator))
        }
        results <- c(sum(results_tmp),sum(results_tmp2)- (sum(results_tmp))^2)
      }
    } else {
      denominator <- sum(c(integrate(function(x) L_binomial(y = 0, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                               lower = -Inf, upper = Inf)$value,
                         integrate(function(x) L_binomial(y = 1, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                   lower = -Inf, upper = Inf)$value))
      numerator <- sum(c(integrate(function(x) x*L_binomial(y = 0, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                             lower = -Inf, upper = Inf)$value,
                         integrate(function(x) x*L_binomial(y = 1, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                                   lower = -Inf, upper = Inf)$value))
      numerator2 <- sum(c(integrate(function(x) x^2*L_binomial(y = 0, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[1],
                      lower = -Inf, upper = Inf)$value,
            integrate(function(x) x^2*L_binomial(y = 1, x, coef = eq[[node]], continuous_part) * prior_gaussian(x, mu_prior, sigma_prior) * p_prior[2],
                      lower = -Inf, upper = Inf)$value))
      results <- c(numerator / denominator, numerator2 / denominator - (numerator / denominator)^2)
    }
  } else {
    lambda_prior <- predictions[[node]]

    eq <- myfit[[children]]
    names(eq) <- sapply(strsplit(names(eq),"[|]"), function(x) x[2])

    bin.nodes <- intersect(names(which(mydists=="binomial")),parents)
    bin.nodes <- setdiff(bin.nodes,node)

    if (length(bin.nodes)>0){
      other.nodes <- parents[-which(parents %in% bin.nodes)]
      other.nodes <- setdiff(other.nodes,node)
    } else {
      other.nodes <- parents
      other.nodes <- setdiff(other.nodes,node)
    }

    predictions_tmp <- predictions[other.nodes]
    predictions_tmp <- lapply(predictions_tmp,function(l){
      l[1]
    })
    continuous_part <- eq[1] + sum(eq[other.nodes]*unlist(predictions_tmp))
    names(continuous_part) <- c()

    if (children %in% names(evidence)){
      p_prior <- c(0,0)
      names(p_prior) <- levels(data[[children]])
      p_prior[grep(predictions[[children]],names(p_prior))] <- 1
    } else {
      p_prior <- predictions[[children]]
    }

    if (length(bin.nodes)>0){
      bin.nodes.evidence <- intersect(names(evidence),bin.nodes)
      if (length(bin.nodes.evidence)>0){
        # at least one bin nodes is an evidence
        predictions_tmp <- predictions[bin.nodes.evidence]
        predictions_tmp <- lapply(predictions_tmp,function(l){
          as.numeric(predictions_tmp[[1]])-1
        })
        continuous_part <- continuous_part + sum(eq[bin.nodes.evidence]*unlist(predictions_tmp))

        if (length(bin.nodes.evidence)==length(bin.nodes)){
          # all bin nodes are evidence
          max_x <- max(1000,4*lambda_prior)
          denominator <- sum(c(sum(sapply(0:max_x, function(x) L_binomial(y = 0, x, coef = eq[[node]], continuous_part) * prior_poisson(x, lambda_prior) * p_prior[1])),
                               sum(sapply(0:max_x, function(x) L_binomial(y = 1, x, coef = eq[[node]], continuous_part) * prior_poisson(x, lambda_prior) * p_prior[2]))))
          numerator <- sum(c(sum(sapply(0:max_x, function(x) x*L_binomial(y = 0, x, coef = eq[[node]], continuous_part) * prior_poisson(x, lambda_prior) * p_prior[1])),
                             sum(sapply(0:max_x, function(x) x*L_binomial(y = 1, x, coef = eq[[node]], continuous_part) * prior_poisson(x, lambda_prior) * p_prior[2]))))

          results <- numerator / denominator
        } else {
          # at least one bin node is not an evidence
          bin.nodes_tmp <- setdiff(bin.nodes,bin.nodes.evidence)
          results_tmp <- c()

          # Generate all combinations of 0 and 1 for binary nodes
          combinations <- expand.grid(rep(list(0:1), length(bin.nodes_tmp)))

          for (i in 1:nrow(combinations)){
            combination_tmp <- as.numeric(combinations[i,])

            # Compute the marginal probability for this combination
            probas <- 1
            for (idx in (1:length(bin.nodes_tmp))) {
              probas <- probas * predictions[[bin.nodes_tmp[idx]]][combination_tmp[idx] + 1]
            }

            # Update the continuous part with contributions from all binary nodes
            continuous_part_tmp <- continuous_part
            for (idx in (1:length(bin.nodes_tmp))) {
              continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes_tmp[idx]] * combination_tmp[idx]
            }

            max_x <- max(1000,4*lambda_prior)
            denominator <- sum(c(sum(sapply(0:max_x, function(x) L_binomial(y = 0, x, coef = eq[[node]], continuous_part_tmp) * prior_poisson(x, lambda_prior) * p_prior[1])),
                                 sum(sapply(0:max_x, function(x) L_binomial(y = 1, x, coef = eq[[node]], continuous_part_tmp) * prior_poisson(x, lambda_prior) * p_prior[2]))))
            numerator <- sum(c(sum(sapply(0:max_x, function(x) x*L_binomial(y = 0, x, coef = eq[[node]], continuous_part_tmp) * prior_poisson(x, lambda_prior) * p_prior[1])),
                               sum(sapply(0:max_x, function(x) x*L_binomial(y = 1, x, coef = eq[[node]], continuous_part_tmp) * prior_poisson(x, lambda_prior) * p_prior[2]))))

            results_tmp <- c(results_tmp,probas*(numerator / denominator))
          }
          results <- sum(results_tmp)
        }
      } else {
        results_tmp <- c()
        # Generate all combinations of 0 and 1 for binary nodes
        combinations <- expand.grid(rep(list(0:1), length(bin.nodes)))

        for (i in 1:nrow(combinations)){
          combination_tmp <- as.numeric(combinations[i,])

          # Compute the marginal probability for this combination
          probas <- 1
          for (idx in (1:length(bin.nodes))) {
            probas <- probas * predictions[[bin.nodes[idx]]][combination_tmp[idx] + 1]
          }

          # Update the continuous part with contributions from all binary nodes
          continuous_part_tmp <- continuous_part
          for (idx in (1:length(bin.nodes))) {
            continuous_part_tmp <- continuous_part_tmp + eq[bin.nodes[idx]] * combination_tmp[idx]
          }

          max_x <- max(1000,4*lambda_prior)
          denominator <- sum(c(sum(sapply(0:max_x, function(x) L_binomial(y = 0, x, coef = eq[[node]], continuous_part_tmp) * prior_poisson(x, lambda_prior) * p_prior[1])),
                               sum(sapply(0:max_x, function(x) L_binomial(y = 1, x, coef = eq[[node]], continuous_part_tmp) * prior_poisson(x, lambda_prior) * p_prior[2]))))
          numerator <- sum(c(sum(sapply(0:max_x, function(x) x*L_binomial(y = 0, x, coef = eq[[node]], continuous_part_tmp) * prior_poisson(x, lambda_prior) * p_prior[1])),
                             sum(sapply(0:max_x, function(x) x*L_binomial(y = 1, x, coef = eq[[node]], continuous_part_tmp) * prior_poisson(x, lambda_prior) * p_prior[2]))))

          results_tmp <- c(results_tmp,probas*(numerator / denominator))
        }
        results <- sum(results_tmp)
      }
    } else {
      max_x <- max(1000,4*lambda_prior)
      denominator <- sum(c(sum(sapply(0:max_x, function(x) L_binomial(y = 0, x, coef = eq[[node]], continuous_part) * prior_poisson(x, lambda_prior) * p_prior[1])),
                           sum(sapply(0:max_x, function(x) L_binomial(y = 1, x, coef = eq[[node]], continuous_part) * prior_poisson(x, lambda_prior) * p_prior[2]))))
      numerator <- sum(c(sum(sapply(0:max_x, function(x) x*L_binomial(y = 0, x, coef = eq[[node]], continuous_part) * prior_poisson(x, lambda_prior) * p_prior[1])),
                         sum(sapply(0:max_x, function(x) x*L_binomial(y = 1, x, coef = eq[[node]], continuous_part) * prior_poisson(x, lambda_prior) * p_prior[2]))))
      results <- numerator / denominator
    }
  }
}

#' Finding the parents of a node in a graph
#'
#' Finds the parents of a given node in a given graph.
#'
#' @param graph A directed graph (igraph object).
#' @param node A node of the graph (either a node label or a number)
#' @return A vector containing the parents of the node (labels of the nodes if the graph is labeled)
#' @examples
#' g <- make_graph("Zachary") # undirected graph
#' find_parents(g, node = 1) # neighbors of the node 1
#' @export
find_parents <- function(graph,node){
  if (is_directed(graph)==FALSE){
    warning(paste0("The provided graph is not directed, this function will output the neighbors of the node ",node,"."))
  }
  if (!is.igraph(graph)){
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

#' Finding the children of a node in a graph
#'
#' Finds the children of a given node in a given graph.
#'
#' @param graph A directed graph (igraph object).
#' @param node A node of the graph (either a node label or a number)
#' @return A vector containing the children of the node (labels of the nodes if the graph is labeled)
#' @examples
#' g <- make_graph("Zachary") # undirected graph
#' find_children(g, node = 1) # neighbors of the node 1
#' @export
find_children <- function(graph,node){
  if (is_directed(graph)==FALSE){
    warning(paste0("The provided graph is not directed, this function will output the neighbors of the node ",node,"."))
  }
  if (!is.igraph(graph)){
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

#' Predicting a root node
#'
#' Predicts a root node in the graph
#'
#' @param data A data frame containing the data (samples in rows, variables in columns).
#' @param mydists A list containing the distributions of the nodes of the graph.
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
predict_root <- function(data, mydists, node){
  if (mydists[[node]]=="gaussian"){
    x <- c(mean(data[[node]]),var(data[[node]]))
  } else if (mydists[[node]]=="poisson"){
    x <- mean(data[[node]])
  } else {
    x <- prop.table(table(data[[node]]))
  }
}

#' Computing the log-likelihood 
#'
#' Computes the log-likelihood of a Poisson variable
#'
#' @param y To fill
#' @param x To fill
#' @param coef To fill
#' @param continuous_part To fill
#' @return To fill
#' @export
LogL_poisson <- function(y, x, coef, continuous_part){
  lambda <- exp(coef*x + continuous_part)
  #dpois(y, lambda)
  #(lambda^y * exp(-lambda)) / gamma(y + 1)
  y * log(lambda) - lambda - lgamma(y + 1)
}

#' Computing the prior
#'
#' Computes the prior of a Gaussian variable
#'
#' @param x To fill
#' @param mu To fill
#' @param sigma2 To fill
#' @return To fill
#' @export
prior_gaussian <- function(x, mu, sigma2){
  dnorm(x, mean=mu, sd=sqrt(sigma2))
}

#' Computing the log-likelihood 
#'
#' Computes the log-likelihood of a Gaussian variable
#'
#' @param y To fill
#' @param x To fill
#' @param coef To fill
#' @param var To fill
#' @param continuous_part To fill
#' @return To fill
#' @export
L_gaussian <- function(y, x, coef, var, continuous_part){
  mu <-  coef*x + continuous_part
  sigma2 <- var
  dnorm(y, mean = mu, sd=sqrt(sigma2))
}

#' Computing the prior
#'
#' Computes the prior of a binomial variable
#'
#' @param x To fill
#' @param p To fill
#' @return To fill
#' @export
prior_binomial <- function(x, p){
  dbinom(x, size=1, prob=p)
}

#' Computing the log-likelihood 
#'
#' Computes the log-likelihood of a binomial variable
#'
#' @param y To fill
#' @param x To fill
#' @param coef To fill
#' @param continuous_part To fill
#' @return To fill
#' @export
L_binomial <- function(y, x, coef, continuous_part){
  mu <- 1 / (1 + exp(-(continuous_part + coef * x )))
  if (y == 1) {
    return(mu)
  } else {
    return(1 - mu)
  }
}

#' Computing the prior
#'
#' Computes the prior of a Poisson variable
#'
#' @param x To fill
#' @param lambda To fill
#' @return To fill
#' @export
prior_poisson <- function(x,lambda){
  dpois(x,lambda)
}

#' Evaluating the performances
#'
#' Evaluates the performances of the predictions
#'
#' @param observations To fill
#' @param predictions To fill
#' @param distribution To fill
#' @param compare.distrib To fill
#' @return To fill
#' @export
EvaluatePerf <- function(observations,predictions,distribution, compare.distrib = FALSE){
  # observations: true values
  # predictions: predicted values
  # distribution: "poisson", "binomial", or "gaussian"
  # plot.distrib: if we want to compare the distributions instead of comparing by samples 
  
  if (compare.distrib == TRUE){
    predictions <- rep(list(predictions), length(observations))
  }
  Scores <- list()
  
  pred.NA <- sapply(predictions, function(l){
    is.na(l[[1]])
  })
  predictions <- predictions[!pred.NA]
  observations <- observations[!pred.NA]
  n <- length(observations)

  if (distribution == "gaussian"){
    # SPE
    SPE <- (sapply(predictions,function(l){l[[1]]})-observations)^2
      
    # log-score
    logScore <- c()
    for (i in (1:length(predictions))){
      logScore_tmp <- 0.5*log(2*pi*predictions[[i]][2]) + (observations[i]-predictions[[i]][1])^2/(2*predictions[[i]][2])
      logScore <- c(logScore,logScore_tmp)
    }
    logScore <- logScore
  
    Scores=c(Scores, list(SPE_mean=mean(SPE),logScore_mean = mean(logScore),SPE = SPE, logScore = logScore))
    
    if (compare.distrib == TRUE){
      # log-likelihood
      logL <- -n/2 * log(2*pi*predictions[[1]][2]) - sum((observations-predictions[[1]][1])^2/(2*predictions[[1]][2]))
      
      # test
      pval <- ks.test(observations,"pnorm", mean=predictions[[1]][1],sd = sqrt(predictions[[1]][2]))$p.value
      
      Scores=c(Scores, list(logL=logL,pval=pval))
    }
  } else if (distribution == "poisson"){
    # SPE
    SPE <- (sapply(predictions,function(l){l[[1]]})-observations)^2
      
    # log-score
    logScore <- c()
    for (i in (1:length(predictions))){
      logScore_tmp <- - observations[i] * log(predictions[[i]]) + predictions[[i]] + lfactorial(observations[[i]])
      logScore <- c(logScore,logScore_tmp)
    }
    logScore <- logScore
    
    Scores=c(Scores, list(SPE_mean=mean(SPE),logScore_mean = mean(logScore),SPE = SPE, logScore = logScore))
      
    if (compare.distrib==TRUE){
      # log-likelihood
      logL <- sum(observations*log(predictions[[1]][1])) - n*log(predictions[[1]][1]) - sum( lfactorial(observations))
        
      # test
      obs_freq <- table(factor(observations, levels = 0:max(observations)))
      expected_freq <- dpois(as.numeric(names(obs_freq)), predictions[[1]][1]) * n
      pval <- chisq.test(x = obs_freq, p = expected_freq / sum(expected_freq))$p.value
      
      Scores=c(Scores, list(logL=logL,pval=pval))
    }
  } else if (distribution == "binomial"){
    # Brier score
    # need to be checked (factor - character,...)
    observations_num <- unlist(sapply(observations, function(l){
      which(l==levels(l)) - 1}))
    BrierScore <- (observations_num-sapply(predictions,function(l){l[[2]]}))^2
  
    # log-score
    # need to be checked
    logScore <- mapply(function(observation,prediction) {
      if (observation == 1) {
        return(log(prediction[2]))  # log(p) if y == 1
      } else {
        return(log(1 - prediction[2]))  # log(1-p) if y == 0
      }
    }, observations,predictions)
  
    logScore <- logScore
      
    Scores=c(Scores, list(BrierScore_mean=mean(BrierScore),logScore_mean = mean(logScore),BrierScore=BrierScore,logScore=logScore))
      
    if (compare.distrib==TRUE){
      # log-likelihood
      logL <- sum(as.numeric(as.character(observations))*log(predictions[[1]][2]) + (1-as.numeric(as.character(observations)))*log(predictions[[1]][1])) 
        
      # test
      observed_counts <- table(observations)
      #expected_counts <- estimation * n
      pval <- chisq.test(observed_counts, p = predictions[[1]], rescale.p = TRUE)$p.value
      
      Scores=c(Scores, list(logL=logL, pval=pval))
    }
  }
  return(Scores)
}

#' Plotting ABN fitted network
#'
#' Plots the ABN network with an emphasis on a specific node
#'
#' @param dag An adjacency matrix (can be the output of the function mostProbable()).
#' @param mydists A list containing the distributions of the nodes of the graph.
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
#' node <- "b4"
#' order <- "up"
#' plot_node <- plot_Abn(dag, mydists, node, order)
#' plot_node
#' @export
plot_Abn <- function (dag, mydists, node, order){ 
  if (!all(colnames(dag) %in% names(mydists))){
    stop("The names of the nodes in dag do not correspond to the ones in mydists.")
  } 
  if (!node %in% colnames(dag)){
    stop("Choose a node that belongs to the network.")
  } 
  if (! order %in% c("up","down")){
    stop("Choose either up or down as an argument for order.")
  }
  mydists <- mydists[colnames(dag)]
  name <- names(mydists)
  
  graph <- graph_from_adjacency_matrix(t(dag))
  
  am.graph <- new(Class = "graphAM", adjMat = dag, edgemode = "directed")
  
  node.shape <- rep(c("circle", "box","ellipse","diamond"), 4)
  shape <- rep(node.shape[1], length(mydists))
  shape[mydists == "binomial"] <- node.shape[2]
  shape[mydists == "poisson"] <- node.shape[3]
  shape[mydists == "multinomial"] <- node.shape[4]
  names(shape) <- names(mydists)
  
  node.fillcolor = c("lightblue", "brown3", "chartreuse3","chartreuse4")
  fillcolor <- rep(node.fillcolor[1], length(mydists))
  names(fillcolor) <- names(mydists)
  if (!is.null(node)) {
    markov.blanket <- abn::mb(dag, node = node, 
                         data.dists = mydists)
    fillcolor[names(mydists) %in% node] <- node.fillcolor[2] # node in red
    
    if (order == "up"){
      parents <- find_parents(graph, node)
      fillcolor[names(mydists) %in% parents] <- node.fillcolor[3]
    } else if (order == "down"){
      children <- find_children(graph, node)
      parents <- unique(unlist(sapply(children,function(l){find_parents(graph, l)})))
      fillcolor[names(mydists) %in% setdiff(parents,node)] <- node.fillcolor[4]
      fillcolor[names(mydists) %in% children] <- node.fillcolor[3]
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

#' Plotting the procedure's workflow
#'
#' Plots and saves an animated gif that represents the procedure's workflow.
#'
#' @param dag An adjacency matrix (can be the output of the function mostProbable()).
#' @param mydists A list containing the distributions of the nodes of the graph.
#' @param hypothesis Node to predict. 
#' @param path A valid path to save the plots
#' @param directory.name The name of the directory that will be created to save the plots.
#' @return An animated gif that represents the procedure's workflow.
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

plot_workflow <- function(dag, mydists, hypothesis, path = NULL, directory.name =NULL){
  path <- createDirectory(path = path, directory.name = directory.name)
  
  graph <- graph_from_adjacency_matrix(t(dag))
  node_order <- names(topo_sort(graph, mode="out"))
  node_max <- which(node_order == hypothesis)
  
  for (i in (1:length(node_order))){
    plots <- plot_Abn(dag, mydists, node = node_order[i], "up")
    png(paste0(path,"/graph",i,".png"), width = 800, height = 600)  
    renderGraph(plots)
    title(main = paste0("Step ",i,", node ",node_order[i],", up"), cex.main = 1.5, font.main = 2)
    dev.off()  
  }
  
  counter <- length(node_order)
  for (i in (length(node_order):node_max)){
    counter <- counter + 1
    plots <- plot_Abn(dag, mydists, node = node_order[i], "down")
    png(paste0(path,"/graph",counter,".png"), width = 800, height = 600) 
    renderGraph(plots)
    title(main = paste0("Step ",counter,", node ",node_order[i],", down"), cex.main = 1.5, font.main = 2)
    dev.off()  
  }
  
  createAnimation(path)
}

#' Creating a directory
#'
#' Creates a directory situated in a particular path 
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

#' Creating an animation 
#'
#' Creates an animated gif file based on png files
#'
#' @param path A valid path to a repository that contains the png to combine in a gif.
#' @return None
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

#' Plotting the posterior distribution
#'
#' Plots the posterior distribution
#'
#' @param predictions To fill.
#' @param hypothesis To fill.
#' @param mydists To fill.
#' @param path To fill.
#' @return None
#' @export
plotPosteriorDistrib <- function(predictions, hypothesis, mydists, path){
  if (mydists[[hypothesis]]=="gaussian"){
    samples <- rnorm(1000, mean = predictions[[hypothesis]][1], sd = sqrt(predictions[[hypothesis]][2]))
    posterior_df <- data.frame(samples)
    
    g <- ggplot(posterior_df, aes(x = samples)) +
      geom_density(fill = "blue", alpha = 0.5) +
      labs(title = paste0("Posterior Distribution of ",hypothesis),
           x = paste0("Value of ",hypothesis),
           y = "Density") +
      theme_minimal()
  } else if (dists[[hypothesis]]=="poisson"){
    values <- 0:(4*(round(predictions[[hypothesis]][1],0)+1))
    posterior_probs <- dpois(values, lambda = predictions[[hypothesis]][1])
    posterior_df <- data.frame(values, posterior_probs)
    
    g <- ggplot(posterior_df, aes(x = values,y=posterior_probs)) +
      geom_bar(stat="identity",fill = "blue", alpha=0.7) +
      scale_x_continuous(breaks = values) +
      labs(title = paste0("Posterior Distribution of ",hypothesis), x = paste0("Value of ",hypothesis), y = "Density") + theme_minimal()
  } else {
    values <- 0:10
    posterior_probs <- dbinom(values, size = length(values), prob = predictions[[hypothesis]][2])
    posterior_df <- data.frame(values, posterior_probs)
    
    g <- ggplot(posterior_df, aes(x = values, y = posterior_probs)) +
      geom_bar(stat = "identity", fill = "blue", alpha = 0.7) +
      scale_x_continuous(breaks = values) +
      labs(title = paste0("Posterior Distribution of ",hypothesis),
           x = "Number of success",
           y = "Probability") +
      theme_minimal()
  }
  ggsave(paste0(path,"/PosteriorDistrib.png"),plot=g,bg="white")
  return(g)
}
