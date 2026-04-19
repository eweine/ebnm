library(testthat)

skip_if_not_installed("ape")

set.seed(1)

tree1 <- ape::read.tree(text = "(t1:1);")
tree2 <- ape::read.tree(text = "(t2:1);")
forest_one_tip <- list(tree1, tree2)

alpha_true <- 1.2
stationary_var_true <- 0.7
mean_var_true <- 0.4

y_one_tip <- c(t1 = 1.5, t2 = -0.8)
s_one_tip <- c(t1 = 0.5, t2 = 0.3)

tree_big <- ape::rtree(18)
tree_big$edge.length <- pmax(tree_big$edge.length, 0.05)
tree_big2 <- ape::rtree(16, tip.label = paste0("clone_tip_", seq_len(16)))
tree_big2$edge.length <- pmax(tree_big2$edge.length, 0.05)
forest_big <- list(tree_big, tree_big2)

simulate_clone_tips <- function(tree, alpha, stationary_var, mean_var) {
  ntip <- ape::Ntip(tree)
  edge_pre <- ape::reorder.phylo(tree, "cladewise")$edge
  elen_pre <- ape::reorder.phylo(tree, "cladewise")$edge.length
  root_node <- setdiff(tree$edge[, 1], tree$edge[, 2])
  mu <- rnorm(1, 0, sqrt(mean_var))
  x_node <- numeric(ntip + tree$Nnode)
  x_node[root_node] <- rnorm(1, mu, sqrt(stationary_var))
  for (i in seq_len(nrow(edge_pre))) {
    p <- edge_pre[i, 1]
    ch <- edge_pre[i, 2]
    a <- exp(-alpha * elen_pre[i])
    q <- stationary_var * (1 - a^2)
    x_node[ch] <- rnorm(1, mu + a * (x_node[p] - mu), sqrt(q))
  }
  setNames(x_node[seq_len(ntip)], tree$tip.label)
}

theta_big <- lapply(
  forest_big,
  simulate_clone_tips,
  alpha = alpha_true,
  stationary_var = stationary_var_true,
  mean_var = mean_var_true
)
s_big <- lapply(theta_big, function(theta) setNames(runif(length(theta), 0.1, 0.3), names(theta)))
y_big <- unlist(Map(function(theta, s) setNames(rnorm(length(theta), theta, s), names(theta)), theta_big, s_big))
s_big_all <- unlist(s_big)

test_that("basic OU clone fit returns a valid ebnm object", {
  fit <- ebnm_tree_ou_clone(y_big, s_big_all, forest_big)
  expect_s3_class(fit, "ebnm")
  expect_s3_class(fit$fitted_g, "tree_ou_clone")
  expect_true(fit$fitted_g$alpha > 0)
  expect_true(fit$fitted_g$stationary_var > 0)
  expect_true(fit$fitted_g$mean_var > 0)
  expect_length(fit$posterior$mean, length(y_big))
  expect_true(is.finite(as.numeric(fit$log_likelihood)))
})

test_that("tree_ou_clone_precomp supports forests", {
  idx <- tree_ou_clone_precomp(forest_big)
  expect_s3_class(idx, "tree_ou_clone_index")
  expect_equal(idx$ntip, length(y_big))
  expect_equal(length(idx$root), 2L)
  expect_length(idx$subindex, 2L)
  expect_length(idx$tree_tip_idx, 2L)
})

test_that("one-tip OU clone forest matches independent normal shrinkage", {
  g_fixed <- tree_ou_clone(alpha_true, stationary_var_true, mean_var_true)
  fit <- ebnm_tree_ou_clone(
    x = y_one_tip,
    s = s_one_tip,
    tree = forest_one_tip,
    g_init = g_fixed,
    fix_g = TRUE
  )

  prior_var <- stationary_var_true + mean_var_true
  expected_var <- stats::setNames(1 / (1 / prior_var + 1 / s_one_tip^2), names(y_one_tip))
  expected_mean <- stats::setNames(expected_var * y_one_tip / s_one_tip^2, names(y_one_tip))

  expect_equal(unname(fit$posterior$mean), unname(expected_mean))
  expect_equal(unname(fit$posterior$sd)^2, unname(expected_var))
  expect_equal(attr(fit$log_likelihood, "df"), 0L)
})

test_that("invalid tree_index is rejected for OU clone fits", {
  expect_error(
    ebnm_tree_ou_clone(y_big, s_big_all, forest_big, tree_index = list(a = 1)),
    "tree_ou_clone_precomp"
  )
})
