library(testthat)

skip_if_not_installed("ape")

set.seed(1)
tree <- ape::rtree(20)
tree$edge.length <- pmax(tree$edge.length, 0.05)
ntip <- ape::Ntip(tree)

alpha_true <- 1.3
stationary_var_true <- 0.7
edge_pre <- ape::reorder.phylo(tree, "cladewise")$edge
elen_pre <- ape::reorder.phylo(tree, "cladewise")$edge.length
root_node <- setdiff(tree$edge[, 1], tree$edge[, 2])
x_node <- numeric(ntip + tree$Nnode)
x_node[root_node] <- rnorm(1, 0, sqrt(stationary_var_true))
for (i in seq_len(nrow(edge_pre))) {
  p <- edge_pre[i, 1]
  ch <- edge_pre[i, 2]
  a <- exp(-alpha_true * elen_pre[i])
  q <- stationary_var_true * (1 - a^2)
  x_node[ch] <- rnorm(1, a * x_node[p], sqrt(q))
}
theta_true <- setNames(x_node[seq_len(ntip)], tree$tip.label)

s_scalar <- 0.2
s_vec <- setNames(runif(ntip, 0.1, 0.3), tree$tip.label)
y <- setNames(rnorm(ntip, theta_true, s_vec), tree$tip.label)

forest <- list(
  tree,
  ape::rtree(12, tip.label = paste0("ou_forest_tip_", seq_len(12)))
)
forest[[2]]$edge.length <- pmax(forest[[2]]$edge.length, 0.05)
forest_s2 <- setNames(runif(12, 0.1, 0.3), forest[[2]]$tip.label)
forest_theta <- lapply(forest, function(tr) {
  tr_pre <- ape::reorder.phylo(tr, "cladewise")
  tr_edge <- tr_pre$edge
  tr_elen <- tr_pre$edge.length
  tr_root <- setdiff(tr$edge[, 1], tr$edge[, 2])
  x_node <- numeric(ape::Ntip(tr) + tr$Nnode)
  x_node[tr_root] <- rnorm(1, 0, sqrt(stationary_var_true))
  for (i in seq_len(nrow(tr_edge))) {
    p <- tr_edge[i, 1]
    ch <- tr_edge[i, 2]
    a <- exp(-alpha_true * tr_elen[i])
    q <- stationary_var_true * (1 - a^2)
    x_node[ch] <- rnorm(1, a * x_node[p], sqrt(q))
  }
  setNames(x_node[seq_len(ape::Ntip(tr))], tr$tip.label)
})
forest_y <- unlist(Map(function(theta, se) {
  setNames(rnorm(length(theta), theta, se), names(theta))
}, forest_theta, list(s_vec, forest_s2)))
forest_s <- c(s_vec, forest_s2)

test_that("basic OU fit returns a valid ebnm object", {
  fit <- ebnm_tree_ou(y, s_vec, tree)
  expect_s3_class(fit, "ebnm")
  expect_s3_class(fit$fitted_g, "tree_ou")
  expect_true(fit$fitted_g$alpha > 0)
  expect_true(fit$fitted_g$stationary_var > 0)
  expect_length(fit$posterior$mean, ntip)
  expect_length(fit$posterior$sd, ntip)
  expect_true(all(fit$posterior$sd > 0))
  expect_true(is.finite(as.numeric(fit$log_likelihood)))
})

test_that("scalar s is expanded correctly for OU fits", {
  fit_scalar <- ebnm_tree_ou(y, s_scalar, tree)
  s_named_uniform <- setNames(rep(s_scalar, ntip), tree$tip.label)
  fit_vec <- ebnm_tree_ou(y, s_named_uniform, tree)
  expect_equal(fit_scalar$fitted_g$alpha, fit_vec$fitted_g$alpha)
  expect_equal(
    fit_scalar$fitted_g$stationary_var,
    fit_vec$fitted_g$stationary_var
  )
  expect_equal(fit_scalar$posterior$mean, fit_vec$posterior$mean)
})

test_that("tree_ou_precomp produces the right class and fields", {
  idx <- tree_ou_precomp(tree)
  expect_s3_class(idx, "tree_ou_index")
  expect_equal(idx$ntip, ntip)
  expect_equal(idx$tip.label, tree$tip.label)
  expect_length(idx$edge_len_to_child, idx$total_nodes)
  expect_length(idx$children, idx$total_nodes)
})

test_that("tree_ou_precomp supports forests", {
  idx <- tree_ou_precomp(forest)
  expect_s3_class(idx, "tree_ou_index")
  expect_equal(idx$ntip, length(forest_y))
  expect_equal(length(idx$root), 2L)
  expect_equal(idx$tip.label, c(tree$tip.label, forest[[2]]$tip.label))
})

test_that("pre-computed OU index gives identical results", {
  idx <- tree_ou_precomp(tree)
  fit_precomp <- ebnm_tree_ou(y, s_vec, tree, tree_index = idx)
  fit_fresh <- ebnm_tree_ou(y, s_vec, tree)
  expect_equal(fit_precomp$fitted_g$alpha, fit_fresh$fitted_g$alpha)
  expect_equal(
    fit_precomp$fitted_g$stationary_var,
    fit_fresh$fitted_g$stationary_var
  )
  expect_equal(fit_precomp$posterior$mean, fit_fresh$posterior$mean)
  expect_equal(fit_precomp$posterior$sd, fit_fresh$posterior$sd)
  expect_equal(
    as.numeric(fit_precomp$log_likelihood),
    as.numeric(fit_fresh$log_likelihood)
  )
})

test_that("fix_g uses supplied OU parameters without optimizing", {
  g_fixed <- tree_ou(alpha_true, stationary_var_true)
  fit_fixed <- ebnm_tree_ou(y, s_vec, tree, g_init = g_fixed, fix_g = TRUE)
  expect_equal(fit_fixed$fitted_g$alpha, alpha_true)
  expect_equal(fit_fixed$fitted_g$stationary_var, stationary_var_true)
  expect_equal(attr(fit_fixed$log_likelihood, "df"), 0L)
})

test_that("stationary root gives the correct one-tip posterior", {
  one_tip_tree <- ape::read.tree(text = "(t1:2);")
  g_fixed <- tree_ou(alpha = 3, stationary_var = 0.8)
  fit <- ebnm_tree_ou(
    x = c(t1 = 1.4),
    s = c(t1 = 0.5),
    tree = one_tip_tree,
    g_init = g_fixed,
    fix_g = TRUE
  )

  expected_var <- 1 / (1 / 0.8 + 1 / 0.25)
  expected_mean <- expected_var * 1.4 / 0.25

  expect_equal(unname(fit$posterior$mean), expected_mean)
  expect_equal(unname(fit$posterior$sd)^2, expected_var)
})

test_that("output is in tree$tip.label order for OU fits", {
  y_shuffled <- y[sample(ntip)]
  s_shuffled <- s_vec[sample(ntip)]
  fit_orig <- ebnm_tree_ou(y, s_vec, tree)
  fit_shuffle <- ebnm_tree_ou(y_shuffled, s_shuffled, tree)
  expect_equal(fit_orig$posterior$mean, fit_shuffle$posterior$mean)
})

test_that("invalid tree_index is rejected for OU fits", {
  expect_error(ebnm_tree_ou(y, s_vec, tree, tree_index = list(a = 1)), "tree_ou_precomp")
})

test_that("forest OU fit matches separate tree fits when g is fixed", {
  g_fixed <- tree_ou(alpha_true, stationary_var_true)
  fit_forest <- ebnm_tree_ou(
    forest_y,
    forest_s,
    forest,
    g_init = g_fixed,
    fix_g = TRUE
  )
  fit_1 <- ebnm_tree_ou(
    forest_y[forest[[1]]$tip.label],
    forest_s[forest[[1]]$tip.label],
    forest[[1]],
    g_init = g_fixed,
    fix_g = TRUE
  )
  fit_2 <- ebnm_tree_ou(
    forest_y[forest[[2]]$tip.label],
    forest_s[forest[[2]]$tip.label],
    forest[[2]],
    g_init = g_fixed,
    fix_g = TRUE
  )

  expect_equal(
    fit_forest$posterior$mean,
    c(fit_1$posterior$mean, fit_2$posterior$mean)
  )
  expect_equal(
    fit_forest$posterior$sd,
    c(fit_1$posterior$sd, fit_2$posterior$sd)
  )
  expect_equal(
    as.numeric(fit_forest$log_likelihood),
    as.numeric(fit_1$log_likelihood) + as.numeric(fit_2$log_likelihood)
  )
})
