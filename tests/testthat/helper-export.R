create_test_abnfit_mle <- function() {
  suppressWarnings({
    # load data
    # load(file = "tests/testthat/testdata/n_250_k_2_groups_Even_mp_7_nodedists_Balanced_s_04_graph1.Rdata")
    load(file = "testdata/n_250_k_2_groups_Even_mp_7_nodedists_Balanced_s_04_graph1.Rdata")

    mycache <- buildScoreCache(data.df = data$data,
                               data.dists = data$dists,
                               method = "mle",
                               max.parents = data$max.parents,
                               group.var = NULL,
                               debugging = FALSE,
                               verbose = FALSE)

    # Structure Learning
    most_prob_dag <- mostProbable(score.cache = mycache, verbose = FALSE)

    # Fit model (this would be actual abn fit)
    fitAbn(object = most_prob_dag, method = "mle")
  })
}

create_test_abnfit_mle_groups <- function() {
  suppressMessages({
    suppressWarnings({
      # load data
      # load(file = "tests/testthat/testdata/n_250_k_2_groups_Even_p_6_nodedists_Balanced_s_04_graph1.rda")
      load(file = "testdata/n_250_k_2_groups_Even_p_6_nodedists_Balanced_s_04_graph1.rda") # no poisson because of fitting issues
      data <- n_250_k_2_groups_Even_p_6_nodedists_Balanced_s_04_graph1$data
      dists <- n_250_k_2_groups_Even_p_6_nodedists_Balanced_s_04_graph1$dists
      max.parents <- n_250_k_2_groups_Even_p_6_nodedists_Balanced_s_04_graph1$max.parents
      groups_samples <- n_250_k_2_groups_Even_p_6_nodedists_Balanced_s_04_graph1$groups_samples

      # Add group variable to data and data dists
      data$groups_samples <- as.factor(groups_samples)

      mycache <- buildScoreCache(data.df = data,
                                 data.dists = dists,
                                 method = "mle",
                                 max.parents = max.parents,
                                 group.var = "groups_samples",
                                 debugging = FALSE,
                                 verbose = FALSE)

      # Structure Learning
      most_prob_dag <- mostProbable(score.cache = mycache, verbose = FALSE)

      # Fit model
      abn_fit <- fitAbn(object = most_prob_dag,
                        method = "mle",
                        group.var = "groups_samples",
                        verbose = FALSE,
                        debugging = FALSE)

      # Store the fitted object for testing because fitting with groups can be slow
      save(abn_fit, file = "tests/testthat/testdata/abnfit_mle_groups.Rdata")
    })
  })
}
