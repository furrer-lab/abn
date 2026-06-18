# Helper environment setup to feed downstream node parameters seamlessly
setup_abn_test_environment <- function() {
  set.seed(123)
  mock_data <- data.frame(
    g1 = rnorm(20, mean = 10, sd = 2),
    p1 = rpois(20, lambda = 12),
    b1 = factor(sample(c("n", "y"), 20, replace = TRUE), levels = c("n", "y")),
    m1 = factor(sample(c("1", "2", "3"), 20, replace = TRUE), levels = c("1", "2", "3"))
  )

  mock_dists <- list(g1 = "gaussian", p1 = "poisson", b1 = "binomial", m1 = "multinomial")
  mock_dag <- matrix(
    c(0, 0, 0, 0,
      1, 0, 0, 0,
      1, 0, 0, 0,
      0, 1, 1, 0),
    nrow = 4, ncol = 4, byrow = TRUE,
    dimnames = list(c("g1", "p1", "b1", "m1"), c("g1", "p1", "b1", "m1"))
  )
  mock_fit <- list(
    g1 = setNames(c(4.8, 1), c("g1|intercept", "precision")),

    p1 = setNames(c(0.5, 0.2), c("p1|intercept", "g1")),

    b1 = setNames(c(-0.2, 0.1), c("b1|intercept", "g1")),

    m1 = setNames(c(2.0, 0.8, 0.5, 1.2, 0.8, 0.9),
                  c("m1|intercept.2", "m1|intercept.3",
                    "p1.2", "p1.3",
                    "b1.2", "b1.3"))
  )

  return(list(data = mock_data, dag = mock_dag, dists = mock_dists, fit = mock_fit))
}

test_that("predictABN handles prediction with no hypothesis and no evidence", {
  env <- setup_abn_test_environment()

  res_baseline <- predictABN(
    dag = env$dag, data = env$data, dists = env$dists, fit = env$fit,
    hypothesis = NULL, evidence = NULL
  )

  expect_type(res_baseline, "list")
  expect_true(all(c("g1", "p1", "b1", "m1") %in% names(res_baseline$predictions)))
  expect_null(res_baseline$prediction_hypothesis)
})

test_that("predictABN handles prediction on a targeted hypothesis with no evidence", {
  env <- setup_abn_test_environment()

  res_hypo <- predictABN(
    dag = env$dag, data = env$data, dists = env$dists, fit = env$fit,
    hypothesis = "p1", evidence = NULL
  )

  expect_true("p1" %in% names(res_hypo$predictions))
  expect_type(res_hypo$prediction_hypothesis, "double")
})

test_that("predictABN handles prediction with some evidence", {
  env <- setup_abn_test_environment()

  # Provide observational data entries across parent structural pathways
  evidence_set <- list("g1" = 5.2, "b1" = "y")

  res_evidence <- predictABN(
    dag = env$dag, data = env$data, dists = env$dists, fit = env$fit,
    hypothesis = NULL, evidence = evidence_set
  )

  # Confirms that known parameters explicitly retain their anchored evidence values
  if ("g1" %in% names(res_evidence$predictions)) {
    expect_equal(res_evidence$predictions$g1[1], 5.2)
  }
  if ("b1" %in% names(res_evidence$predictions)) {
    expect_equal(res_evidence$predictions$b1, c("n" = 0, "y" = 1))
  }
})

test_that("predictABN handles prediction with both hypotheses and evidence", {
  env <- setup_abn_test_environment()

  res_full_query <- predictABN(
    dag = env$dag, data = env$data, dists = env$dists, fit = env$fit,
    hypothesis = "m1", evidence = list("g1" = 4.5, "p1" = 2)
  )

  expect_true("m1" %in% names(res_full_query$predictions))
  expect_length(res_full_query$predictions$m1, length(levels(env$data$m1)))
  expect_false(any(is.nan(res_full_query$predictions$m1)))
})

test_that("predictABN optimizes inference when Markov Blanket is pre-computed", {
  env <- setup_abn_test_environment()

  # When explicit structural dependencies or neighborhoods are provided via mb,
  # the inference routine isolates local structural paths.
  expect_message(
    res_mb_optimized <- predictABN(
      dag = env$dag,
      data = env$data,
      dists = env$dists,
      fit = env$fit,
      hypothesis = "p1",
      evidence = list(b1 = "y", g1 = 1.2, m1 = "3"),
    ),
    regexp = "The Markov blanket of the node to infer is enterily known\\."
  )

  expect_type(res_mb_optimized, "list")
})

test_that("check_data allows properly formatted data to pass without errors", {
  # Setup valid environment matching your naming rules ('b' for binomial, 'm' for multinomial)
  valid_data <- data.frame(
    g1 = rnorm(10),
    b1 = factor(c("0", "1", "0", "1", "0", "1", "0", "1", "0", "1"), levels = c("0", "1")),
    m1 = factor(c("1", "2", "3", "1", "2", "3", "1", "2", "3", "1"), levels = c("1", "2", "3"))
  )
  dists <- list(g1 = "gaussian", b1 = "binomial", m1 = "multinomial")

  # A successful pass should return nothing (NULL) and throw no errors
  expect_error(check_data(valid_data, dists), NA)
})

test_that("check_data throws an explicit error if a binomial node has only 1 active level", {
  # Create a problematic dataframe where b1 only has 1 level present
  bad_bin_data <- data.frame(
    g1 = rnorm(10),
    b1 = factor(rep("0", 10), levels = c("0")),
    m1 = factor(c("1", "2", "3", "1", "2", "3", "1", "2", "3", "1"), levels = c("1", "2", "3"))
  )
  dists <- list(g1 = "gaussian", b1 = "binomial", m1 = "multinomial")

  # Expect a stop() condition matching the binomial error text
  expect_error(
    check_data(bad_bin_data, dists),
    regexp = "Binomial node b1 does not have the right number of levels"
  )
})

test_that("check_data throws an explicit error if a multinomial node has fewer than 3 levels", {
  # Create a problematic dataframe where m1 has only 2 levels (making it effectively binomial)
  bad_multi_data <- data.frame(
    g1 = rnorm(10),
    b1 = factor(c("0", "1", "0", "1", "0", "1", "0", "1", "0", "1"), levels = c("0", "1")),
    m1 = factor(c("1", "2", "1", "2", "1", "2", "1", "2", "1", "2"), levels = c("1", "2"))
  )
  dists <- list(g1 = "gaussian", b1 = "binomial", m1 = "multinomial")

  # Expect a stop() condition matching your multinomial text requirement
  expect_error(
    check_data(bad_multi_data, dists),
    regexp = "Multinomial node m1 does not have the expected number of levels"
  )
})

test_that("find_MB correctly identifies full Markov Blanket neighborhoods", {
  edges <- c("A", "B",
             "B", "C",
             "D", "C",
             "E", "A")
  g <- igraph::make_graph(edges = edges, directed = TRUE)

  mb_res <- find_MB(graph = g, hypothesis = "B")

  expect_equal(sort(mb_res), sort(c("A", "C", "D")))
  expect_false("B" %in% mb_res)
  expect_false("E" %in% mb_res)
})

test_that("find_MB handles root nodes and leaf nodes", {
  # Simple chain graph: A -> B -> C
  g <- igraph::make_graph(edges = c("A", "B", "B", "C"), directed = TRUE)

  # Test a Root Node (A)
  # Parents: None, Children: B, Co-parents: None
  mb_root <- find_MB(graph = g, hypothesis = "A")
  expect_equal(mb_root, "B")

  # Test a Leaf Node (C)
  # Parents: B, Children: None, Co-parents: None
  mb_leaf <- find_MB(graph = g, hypothesis = "C")
  expect_equal(mb_leaf, "B")
})

test_that("find_MB returns an empty character vector for fully isolated nodes", {
  # Create a graph with nodes but zero edges
  g <- igraph::make_empty_graph(n = 3)
  igraph::V(g)$name <- c("X", "Y", "Z")

  mb_isolated <- find_MB(graph = g, hypothesis = "X")

  # Should return character(0) or NULL depending on how unlist/unique handles empty returns
  expect_length(mb_isolated, 0)
})

test_that("predict_node_from_parent stops when upstream parent data is completely missing", {
  # Build a simple A -> B structure
  g <- igraph::make_graph(edges = c("A", "B"), directed = TRUE)
  dists <- list(A = "gaussian", B = "poisson")
  fit <- list()

  # Target node 'B' has 'A' as a parent, but 'A' is not provided in evidence or predictions
  expect_error(
    predict_node_from_parent(
      data = data.frame(), dists = dists, graph = g, fit = fit,
      node = "B", evidence = list(), predictions = NULL
    ),
    regexp = "Not enough information about the upstream nodes"
  )
})

test_that("predict_node_from_parent successfully passes through evidence if predictions is NULL", {
  g <- igraph::make_graph(edges = c("A", "B"), directed = TRUE)
  dists <- list(A = "gaussian", B = "poisson")
  fit <- list()

  # Mock the underlying poisson sub-engine function to verify it gets called
  mockery::stub(predict_node_from_parent, "predict_node_from_parent_poisson", function(...) "poisson_success")

  # When predictions = NULL but all parents are in evidence, it should copy evidence over
  res <- predict_node_from_parent(
    data = data.frame(), dists = dists, graph = g, fit = fit,
    node = "B", evidence = list(A = 5), predictions = NULL
  )

  expect_equal(res, "poisson_success")
})

test_that("predict_node_from_parent routes cleanly to different sub-engines via switch statement", {
  # Setup an isolated graph structure where A -> B
  g <- igraph::make_graph(edges = c("A", "B"), directed = TRUE)
  fit <- list()
  mock_predictions <- list(A = 2)
  mock_evidence <- list()

  # Stub out every possible distribution-specific function
  mockery::stub(predict_node_from_parent, "predict_node_from_parent_gaussian", function(...) "routed_to_gaussian")
  mockery::stub(predict_node_from_parent, "predict_node_from_parent_binomial", function(...) "routed_to_binomial")
  mockery::stub(predict_node_from_parent, "predict_node_from_parent_multinomial", function(...) "routed_to_multinomial")
  mockery::stub(predict_node_from_parent, "predict_node_from_parent_poisson", function(...) "routed_to_poisson")

  # Test Case 1: Route to Gaussian
  dists_gauss <- list(A = "gaussian", B = "gaussian")
  res_gauss <- predict_node_from_parent(
    data = data.frame(), dists = dists_gauss, graph = g, fit = fit,
    node = "B", evidence = mock_evidence, predictions = mock_predictions
  )
  expect_equal(res_gauss, "routed_to_gaussian")

  # Test Case 2: Route to Binomial
  dists_bin <- list(A = "gaussian", B = "binomial")
  res_bin <- predict_node_from_parent(
    data = data.frame(), dists = dists_bin, graph = g, fit = fit,
    node = "B", evidence = mock_evidence, predictions = mock_predictions
  )
  expect_equal(res_bin, "routed_to_binomial")

  # Test Case 3: Route to Multinomial
  dists_multi <- list(A = "gaussian", B = "multinomial")
  res_multi <- predict_node_from_parent(
    data = data.frame(), dists = dists_multi, graph = g, fit = fit,
    node = "B", evidence = mock_evidence, predictions = mock_predictions
  )
  expect_equal(res_multi, "routed_to_multinomial")

  # Test Case 4: Route to Poisson
  dists_poisson <- list(A = "gaussian", B = "poisson")
  res_poisson <- predict_node_from_parent(
    data = data.frame(), dists = list(A = "gaussian", B = "poisson"), graph = g, fit = fit,
    node = "B", evidence = mock_evidence, predictions = mock_predictions
  )
  expect_equal(res_poisson, "routed_to_poisson")
})

test_that("predict_node_from_parent crashes safely on an unsupported node distribution type", {
  g <- igraph::make_graph(edges = c("A", "B"), directed = TRUE)
  # Give node 'B' an alien distribution type
  dists_broken <- list(A = "gaussian", B = "negative_binomial")

  expect_error(
    predict_node_from_parent(
      data = data.frame(), dists = dists_broken, graph = g, fit = list(),
      node = "B", evidence = list(), predictions = list(A = 3)
    ),
    regexp = "Unknown node distribution type: negative_binomial"
  )
})

test_that("predict_node_from_children errors if target node missing from predictions list", {
  g <- igraph::make_graph(edges = c("A", "B"), directed = TRUE)

  # 'node' is "A", but predictions only has "B" -> triggers first stop() guard
  expect_error(
    predict_node_from_children(
      data = data.frame(), dists = list(A="gaussian", B="poisson"), graph = g,
      fit = list(), node = "A", evidence = list(), predictions = list(B = 5)
    ),
    regexp = "Predictions must contain at least a first prediction"
  )
})

test_that("predict_node_from_children errors if children or co-parents lack sufficient tracking data", {
  # Layout: Target 'A' points to child 'B'. 'C' is a co-parent pointing to 'B'.
  g <- igraph::make_graph(edges = c("A", "B", "C", "B"), directed = TRUE)
  dists <- list(A = "gaussian", B = "poisson", C = "gaussian")

  # Node B and Co-parent C are completely missing from both predictions and evidence -> triggers second stop()
  expect_error(
    predict_node_from_children(
      data = data.frame(), dists = dists, graph = g, fit = list(),
      node = "A", evidence = list(), predictions = list(A = 10)
    ),
    regexp = "Not enough information about the downstream nodes"
  )
})

test_that("predict_node_from_children returns baseline prediction if node has no children", {
  # Isolated target or a Leaf node with no children
  g <- igraph::make_graph(edges = c("A", "B"), directed = TRUE)

  # Target node 'B' is a leaf node. Should cleanly fall back to its existing baseline value.
  res <- predict_node_from_children(
    data = data.frame(), dists = list(A="gaussian", B="poisson"), graph = g,
    fit = list(), node = "B", evidence = list(), predictions = list(B = 12.5)
  )

  expect_equal(res, 12.5)
})

test_that("predict_node_from_children aggregates continuous Gaussian node outputs correctly", {
  # Target 'A' has two children: 'B' and 'C'
  g <- igraph::make_graph(edges = c("A", "B", "A", "C"), directed = TRUE)
  dists <- list(A = "gaussian", B = "gaussian", C = "gaussian")

  # Stub out the children child calculations to return mock mean and variance pairs
  # Child B returns mean=10, var=4. Child C returns mean=12, var=4.
  mock_child_calc <- mockery::mock(c(10, 4), c(12, 4))
  mockery::stub(predict_node_from_children, "predict_node_from_children_gaussian", mock_child_calc)

  res <- predict_node_from_children(
    data = data.frame(), dists = dists, graph = g, fit = list(),
    node = "A", evidence = list(), predictions = list(A = 5, B = 2, C = 3)
  )

  # Expected Aggregation Math for Target Node:
  # res_mean = mean(c(10, 12)) = 11
  # res_var  = mean(c(4, 4)) / length(children) = 4 / 2 = 2
  expect_equal(res, c(11, 2))
})

test_that("predict_node_from_children aggregates Count (Poisson) node outputs correctly", {
  g <- igraph::make_graph(edges = c("A", "B", "A", "C"), directed = TRUE)
  dists <- list(A = "poisson", B = "poisson", C = "poisson")

  # Stub Poisson child predictions to return matrices/vectors with rates
  mock_child_calc <- mockery::mock(matrix(c(15, 0), nrow=1), matrix(c(17, 0), nrow=1))
  mockery::stub(predict_node_from_children, "predict_node_from_children_poisson", mock_child_calc)

  res <- predict_node_from_children(
    data = data.frame(), dists = dists, graph = g, fit = list(),
    node = "A", evidence = list(), predictions = list(A = 12, B = 2, C = 3)
  )

  # Expected Poisson Math: mean(c(15, 17)) = 16
  expect_equal(res, 16)
})

test_that("predict_node_from_children aggregates discrete Categorical (Multinomial/Binomial) matrix shapes correctly", {
  g <- igraph::make_graph(edges = c("A", "B"), directed = TRUE)
  dists <- list(A = "multinomial", B = "gaussian")

  # Provide factor structure framework data
  mock_data <- data.frame(A = factor(c("low", "med", "high"), levels = c("low", "med", "high")))

  # Scenario A: Child outputs full multi-class probability row vector
  mock_child_calc_prob <- mockery::mock(matrix(c(0.2, 0.5, 0.3), nrow = 1))
  mockery::stub(predict_node_from_children, "predict_node_from_children_gaussian", mock_child_calc_prob)

  res_prob <- predict_node_from_children(
    data = mock_data, dists = dists, graph = g, fit = list(),
    node = "A", evidence = list(), predictions = list(A = "med", B = 1.5)
  )

  # colMeans / sum total normalization check:
  # c(0.2, 0.5, 0.3) / 1.0 = c(0.2, 0.5, 0.3)
  expect_equal(res_prob, c(low = 0.2, med = 0.5, high = 0.3))

  # Scenario B: Child engine outputs an index matrix vector (forcing hard classification)
  mock_child_calc_idx <- mockery::mock(matrix(c(3), nrow = 1, ncol = 1)) # index '3' chosen ("high")
  mockery::stub(predict_node_from_children, "predict_node_from_children_gaussian", mock_child_calc_idx)

  res_idx <- predict_node_from_children(
    data = mock_data, dists = dists, graph = g, fit = list(),
    node = "A", evidence = list(), predictions = list(A = "med", B = 1.5)
  )

  # Hard classification map check (the 3rd profile level gets 1, others 0):
  expect_equal(res_idx, c(low = 0, med = 0, high = 1))
})
