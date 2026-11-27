create_test_abnfit_mle <- function() {
  load("testdata/n_250_k_2_groups_Even_mp_7_nodedists_Balanced_s_04_graph1.Rdata")

  suppressWarnings({
    dag <- matrix(
      0,
      nrow = length(data$dists), ncol = length(data$dists),
      dimnames = list(names(data$dists), names(data$dists))
    )

    # Fixed, deliberately simple DAG:
    #   g2 <- g1, p2 <- p1, m1 and m2 are roots.
    # No multinomial parents and no parent-set search, so the design is stable
    # across numerical environments.
    dag["g2", "g1"] <- 1
    dag["p2", "p1"] <- 1

    fitAbn(
      dag = dag,
      data.df = data$data,
      data.dists = data$dists,
      method = "mle"
    )
  })
}

# `testdata/abnfit_mle_groups.Rdata` was generated from
# `testdata/n_250_k_2_groups_Even_p_6_nodedists_Balanced_s_04_graph1.rda`
# with `buildScoreCache(..., group.var = "groups_samples")`, `mostProbable()`,
# and `fitAbn(..., group.var = "groups_samples")`.
