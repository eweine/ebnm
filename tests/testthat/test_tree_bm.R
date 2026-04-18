library(testthat)

skip_if_not_installed("ape")

set.seed(1)
tree <- ape::rtree(20)
tree$edge.length <- pmax(tree$edge.length, 0.05)
ntip <- ape::Ntip(tree)

# Simulate tip states under BM with root fixed at 0
bm_var_true <- 0.5
edge_pre <- ape::reorder.phylo(tree, "cladewise")$edge
elen_pre  <- ape::reorder.phylo(tree, "cladewise")$edge.length
root_node <- setdiff(tree$edge[, 1], tree$edge[, 2])
x_node    <- rep(0, ntip + tree$Nnode)
for (i in seq_len(nrow(edge_pre))) {
  p <- edge_pre[i, 1]; ch <- edge_pre[i, 2]
  x_node[ch] <- rnorm(1, x_node[p], sqrt(bm_var_true * elen_pre[i]))
}
theta_true <- setNames(x_node[seq_len(ntip)], tree$tip.label)

s_scalar <- 0.2
s_vec    <- setNames(runif(ntip, 0.1, 0.3), tree$tip.label)
y        <- setNames(rnorm(ntip, theta_true, s_vec), tree$tip.label)

test_that("basic fit returns a valid ebnm object", {
  fit <- ebnm_tree_bm(y, s_vec, tree)
  expect_s3_class(fit, "ebnm")
  expect_s3_class(fit$fitted_g, "tree_bm")
  expect_true(fit$fitted_g$bm_var > 0)
  expect_length(fit$posterior$mean, ntip)
  expect_length(fit$posterior$sd, ntip)
  expect_true(all(fit$posterior$sd > 0))
  expect_true(is.finite(as.numeric(fit$log_likelihood)))
})

test_that("scalar s is expanded correctly", {
  fit_scalar <- ebnm_tree_bm(y, s_scalar, tree)
  s_named_uniform <- setNames(rep(s_scalar, ntip), tree$tip.label)
  fit_vec    <- ebnm_tree_bm(y, s_named_uniform, tree)
  expect_equal(fit_scalar$fitted_g$bm_var, fit_vec$fitted_g$bm_var)
  expect_equal(fit_scalar$posterior$mean, fit_vec$posterior$mean)
})

test_that("s vector of wrong length errors", {
  expect_error(ebnm_tree_bm(y, c(0.1, 0.2), tree), "length")
})

test_that("tree_bm_precomp produces the right class and fields", {
  idx <- tree_bm_precomp(tree)
  expect_s3_class(idx, "tree_bm_index")
  expect_equal(idx$ntip, ntip)
  expect_equal(idx$tip.label, tree$tip.label)
  expect_length(idx$edge_len_to_child, idx$total_nodes)
  expect_length(idx$children, idx$total_nodes)
})

test_that("pre-computed index gives identical results to building on the fly", {
  idx <- tree_bm_precomp(tree)
  fit_precomp <- ebnm_tree_bm(y, s_vec, tree, tree_index = idx)
  fit_fresh   <- ebnm_tree_bm(y, s_vec, tree)
  expect_equal(fit_precomp$fitted_g$bm_var, fit_fresh$fitted_g$bm_var)
  expect_equal(fit_precomp$posterior$mean,  fit_fresh$posterior$mean)
  expect_equal(fit_precomp$posterior$sd,    fit_fresh$posterior$sd)
  expect_equal(as.numeric(fit_precomp$log_likelihood),
               as.numeric(fit_fresh$log_likelihood))
})

test_that("fix_g uses supplied bm_var without optimizing", {
  g_fixed <- tree_bm(bm_var_true)
  fit_fixed <- ebnm_tree_bm(y, s_vec, tree, g_init = g_fixed, fix_g = TRUE)
  expect_equal(fit_fixed$fitted_g$bm_var, bm_var_true)
  expect_equal(attr(fit_fixed$log_likelihood, "df"), 0L)
})

test_that("output argument controls returned fields", {
  fit_min <- ebnm_tree_bm(y, s_vec, tree, output = c("posterior_mean", "fitted_g"))
  expect_null(fit_min$data)
  expect_null(fit_min$log_likelihood)
  expect_false(is.null(fit_min$posterior))
  expect_false(is.null(fit_min$fitted_g))

  fit_all <- suppressWarnings(ebnm_tree_bm(y, s_vec, tree, output = ebnm_output_all()))
  expect_false(is.null(fit_all$data))
  expect_false(is.null(fit_all$posterior$lfsr))
})

test_that("output is in tree$tip.label order", {
  y_shuffled <- y[sample(ntip)]
  s_shuffled <- s_vec[sample(ntip)]
  fit_orig    <- ebnm_tree_bm(y, s_vec, tree)
  fit_shuffle <- ebnm_tree_bm(y_shuffled, s_shuffled, tree)
  expect_equal(fit_orig$posterior$mean, fit_shuffle$posterior$mean)
})

test_that("invalid tree_index is rejected", {
  expect_error(ebnm_tree_bm(y, s_vec, tree, tree_index = list(a = 1)), "tree_bm_precomp")
})
