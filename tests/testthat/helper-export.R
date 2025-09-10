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
