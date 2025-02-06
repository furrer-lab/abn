predictABN <- function(data, mydists, dag, myfit, hypothesis, evidence){
  # main function to predict in ABN
  # inputs:
  # data: the dataset
  # mydists: distributions of the node
  # dag: the DAG (can be the output of function mostProbable mp.dag$dag
  # myfit: output of function fitAbn (myfit$modes)
  # hypothesis: node to predict
  # evidence: known nodes

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
  return(list(prediction_hypothesis = predictions[[hypothesis]], predictions = predictions))
}

check_evidence <- function(data, mydists, hypothesis, evidence){
  # this function check the format of the evidence and eventually clean it
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

predict_node_from_parent_poisson <- function(data, mydists, myfit, node, evidence, parents, predictions){
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

predict_node_from_parent_gaussian <- function(data, mydists, myfit, node, evidence, parents, predictions){
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

predict_node_from_parent_binomial <- function(data, mydists, myfit, node, evidence, parents, predictions){
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
          as.numeric(levels(predictions_tmp[[1]])[as.numeric(predictions_tmp[[1]])])
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
          as.numeric(levels(predictions_tmp[[1]])[as.numeric(predictions_tmp[[1]])])
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
            results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator)-(probas*(numerator / denominator))^2)
          }
          results <- c(sum(results_tmp),sum(results_tmp2))
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
          results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator)-(probas*(numerator / denominator))^2)
        }
        results <- c(sum(results_tmp),sum(results_tmp2))
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
          as.numeric(levels(predictions_tmp[[1]])[as.numeric(predictions_tmp[[1]])])
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
          as.numeric(levels(predictions_tmp[[1]])[as.numeric(predictions_tmp[[1]])])
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
            } else {
              denominator <- try.denominator$value
            }

            try.numerator <- try(integrate(function(x) x*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                           lower = -Inf, upper = Inf),TRUE)
            if (length(try.numerator)==1){
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
            } else {
              numerator2 <- try.numerator2$value
            }
            results_tmp <- c(results_tmp,probas*(numerator / denominator))
            results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator)-(probas*(numerator / denominator))^2)
          }
          results <- c(sum(results_tmp),sum(results_tmp2))
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
            numerator2 <- integrate(function(x) x^2*exp(LogL_poisson(y = predictions[[children]], x, coef = eq[[node]], continuous_part_tmp)) * prior_gaussian(x, mu_prior, sigma_prior),
                                   lower = -10, upper = 10)$value
          } else {
            numerator2 <- try.numerator2$value
          }
          
          results_tmp <- c(results_tmp,probas*(numerator / denominator))
          results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator)-(probas*(numerator / denominator))^2)
        }
        results <- c(sum(results_tmp),sum(results_tmp2))
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
          as.numeric(levels(predictions_tmp[[1]])[as.numeric(predictions_tmp[[1]])])
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
          as.numeric(levels(predictions_tmp[[1]])[as.numeric(predictions_tmp[[1]])])
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

    if (length(bin.nodes)>0){
      bin.nodes.evidence <- intersect(names(evidence),bin.nodes)

      if (length(bin.nodes.evidence)>0){
        # at least one bin nodes is an evidence
        predictions_tmp <- predictions[bin.nodes.evidence]
        predictions_tmp <- lapply(predictions_tmp,function(l){
          as.numeric(levels(predictions_tmp[[1]])[as.numeric(predictions_tmp[[1]])])
        })
        continuous_part <- continuous_part + sum(eq[bin.nodes.evidence]*unlist(predictions_tmp))

        if (length(bin.nodes.evidence)==length(bin.nodes)){
          # all bin nodes are evidence
          numerator <-  function(x,y){
            L_binomial(y, x, coef = eq[[node]], continuous_part)  * prior_binomial(x, p_prior[2])
          }

          numerator_2 <- function(x){
            numerator(0,x)*p_prior[1] + numerator(1,x)*p_prior[2]
          }

          denominator <- numerator_2(0) + numerator_2(1)
          results <- c(numerator_2(0) / denominator,numerator_2(1) / denominator)
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

            numerator_2 <- function(x){
              numerator(0,x)*p_prior[1] + numerator(1,x)*p_prior[2]
            }

            denominator <- numerator_2(0) + numerator_2(1)
            results_tmp <- c(results_tmp,probas*numerator_2(0) / denominator)
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

          numerator_2 <- function(x){
            numerator(0,x)*p_prior[1] + numerator(1,x)*p_prior[2]
          }

          denominator <- numerator_2(0) + numerator_2(1)
          results_tmp <- c(results_tmp,probas*numerator_2(0) / denominator)
        }
        results <- c(sum(results_tmp),1-sum(results_tmp))
      }
    } else {
      numerator <-  function(x,y){
        L_binomial(y, x, coef = eq[[node]], continuous_part)  * prior_binomial(x, p_prior[2])
      }

      numerator_2 <- function(x){
        numerator(0,x)*p_prior[1] + numerator(1,x)*p_prior[2]
      }

      denominator <- numerator_2(0) + numerator_2(1)
      results <- c(numerator_2(0) / denominator,numerator_2(1) / denominator)
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
          as.numeric(levels(predictions_tmp[[1]])[as.numeric(predictions_tmp[[1]])])
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
            results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator)-(probas*(numerator / denominator))^2)
          }
          results <- c(sum(results_tmp),sum(results_tmp2))
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
          results_tmp2 <- c(results_tmp2,probas*(numerator2 / denominator)-(probas*(numerator / denominator))^2)
        }
        results <- c(sum(results_tmp),sum(results_tmp2))
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
          as.numeric(levels(predictions_tmp[[1]])[as.numeric(predictions_tmp[[1]])])
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

find_parents <- function(graph,node){
  # find the parents of the node of interest
  parents <- names(neighbors(graph, v= node,mode="in"))
  return(parents)
}

find_children <- function(graph,node){
  # find the children of the node of interest
  children <- names(neighbors(graph, v= node,mode="out"))
  return(children)
}

predict_root <- function(data, mydists, node){
  if (mydists[[node]]=="gaussian"){
    x <- c(mean(data[[node]]),var(data[[node]]))
  } else if (mydists[[node]]=="poisson"){
    x <- mean(data[[node]])
  } else {
    x <- prop.table(table(data[[node]]))
  }
}

LogL_poisson <- function(y, x, coef, continuous_part){
  lambda <- exp(coef*x + continuous_part)
  #dpois(y, lambda)
  #(lambda^y * exp(-lambda)) / gamma(y + 1)
  y * log(lambda) - lambda - lgamma(y + 1)
}

prior_gaussian <- function(x, mu, sigma2){
  dnorm(x, mean=mu, sd=sqrt(sigma2))
}

L_gaussian <- function(y, x, coef, var, continuous_part){
  mu <-  coef*x + continuous_part
  sigma2 <- var
  dnorm(y, mean = mu, sd=sqrt(sigma2))
}

prior_binomial <- function(x, p){
  dbinom(x, size=1, prob=p)
}

L_binomial <- function(y, x, coef, continuous_part){
  mu <- 1 / (1 + exp(-(continuous_part + coef * x )))
  if (y == 1) {
    return(mu)
  } else {
    return(1 - mu)
  }
}

prior_poisson <- function(x,lambda){
  dpois(x,lambda)
}

EvaluatePerf2 <- function(observations,predictions,distribution){
  # observations: true values
  # predictions: predicted values
  # distribution: "poisson", "binomial", or "gaussian"
  pred.NA <- sapply(predictions, function(l){
    is.na(l[[1]])
  })
  predictions <- predictions[!pred.NA]
  observations <- observations[!pred.NA]
  n <- length(observations)

  if (distribution == "binomial"){
    # first method with AUC

    # second method with Brier score
    BrierScore <- (1/n)*sum((as.numeric(observations)-sapply(predictions,function(l){l[[2]]}))^2)

    # third method with log-score
    LogScore <- mapply(function(observation,prediction) {
      if (observation == 1) {
        return(log(prediction[2]))  # log(p) if y == 1
      } else {
        return(log(1 - prediction[2]))  # log(1-p) if y == 0
      }
    }, observations,predictions)

    LogScore <- -mean(LogScore)

    Scores <- c(BrierScore,LogScore)
    names(Scores) <- c("BrierScore","LogScore")
  } else {
    # first method with squared prediction error
    SPE <- mean((sapply(predictions,function(l){l[[1]]})-observations)^2)

    # second method with log-score
    LogScore <- -mean(sapply(predictions,function(l){log(l[[2]])}))

    Scores <- c(SPE,LogScore)
    names(Scores) <- c("SPE","LogScore")
  }
  return(Scores)
}
