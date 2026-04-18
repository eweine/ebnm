library(ape)

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", helpers = FALSE, quiet = TRUE)
} else {
  library(ebnm)
}

set.seed(1)

make_tree <- function(n_tip, prefix) {
  tr <- ape::rtree(n_tip, tip.label = paste0(prefix, seq_len(n_tip)))
  tr$edge.length <- pmax(tr$edge.length, 0.05)
  tr
}

simulate_bm_tips <- function(tree, bm_var) {
  ntip <- ape::Ntip(tree)
  edge_pre <- ape::reorder.phylo(tree, "cladewise")$edge
  elen_pre <- ape::reorder.phylo(tree, "cladewise")$edge.length
  x_node <- numeric(ntip + tree$Nnode)

  for (i in seq_len(nrow(edge_pre))) {
    p <- edge_pre[i, 1]
    ch <- edge_pre[i, 2]
    x_node[ch] <- rnorm(1, x_node[p], sqrt(bm_var * elen_pre[i]))
  }

  setNames(x_node[seq_len(ntip)], tree$tip.label)
}

simulate_ou_tips <- function(tree, alpha, stationary_var) {
  ntip <- ape::Ntip(tree)
  edge_pre <- ape::reorder.phylo(tree, "cladewise")$edge
  elen_pre <- ape::reorder.phylo(tree, "cladewise")$edge.length
  root_node <- setdiff(tree$edge[, 1], tree$edge[, 2])
  x_node <- numeric(ntip + tree$Nnode)
  x_node[root_node] <- rnorm(1, 0, sqrt(stationary_var))

  for (i in seq_len(nrow(edge_pre))) {
    p <- edge_pre[i, 1]
    ch <- edge_pre[i, 2]
    a <- exp(-alpha * elen_pre[i])
    q <- stationary_var * (1 - a^2)
    x_node[ch] <- rnorm(1, a * x_node[p], sqrt(q))
  }

  setNames(x_node[seq_len(ntip)], tree$tip.label)
}

fit_separately <- function(trees, y_list, s_list, fit_fun, g_fixed) {
  fits <- Map(function(tr, y, s) {
    fit_fun(y, s, tr, g_init = g_fixed, fix_g = TRUE)
  }, trees, y_list, s_list)

  list(
    posterior_mean = unlist(lapply(fits, function(fit) fit$posterior$mean)),
    posterior_sd = unlist(lapply(fits, function(fit) fit$posterior$sd)),
    loglik = sum(vapply(fits, function(fit) as.numeric(fit$log_likelihood), numeric(1)))
  )
}

tree_list <- list(
  make_tree(12, "a_"),
  make_tree(10, "b_"),
  make_tree(14, "c_")
)

cat("forest sizes:", paste(vapply(tree_list, ape::Ntip, integer(1)), collapse = ", "), "\n")

bm_var_true <- 0.7
bm_theta <- lapply(tree_list, simulate_bm_tips, bm_var = bm_var_true)
bm_s <- lapply(bm_theta, function(theta) setNames(runif(length(theta), 0.1, 0.3), names(theta)))
bm_y_list <- Map(function(theta, s) {
  setNames(rnorm(length(theta), theta, s), names(theta))
}, bm_theta, bm_s)
bm_y <- unlist(bm_y_list)
bm_s_all <- unlist(bm_s)

bm_forest_fit <- ebnm_tree_bm(
  bm_y,
  bm_s_all,
  tree_list,
  g_init = tree_bm(bm_var_true),
  fix_g = TRUE
)
bm_sep_fit <- fit_separately(tree_list, bm_y_list, bm_s, ebnm_tree_bm, tree_bm(bm_var_true))

cat("\nBM forest check\n")
cat(
  "max abs diff posterior mean:",
  signif(max(abs(bm_forest_fit$posterior$mean - bm_sep_fit$posterior_mean)), 4),
  "\n"
)
cat(
  "max abs diff posterior sd:",
  signif(max(abs(bm_forest_fit$posterior$sd - bm_sep_fit$posterior_sd)), 4),
  "\n"
)
cat(
  "forest loglik - sum separate:",
  signif(as.numeric(bm_forest_fit$log_likelihood) - bm_sep_fit$loglik, 4),
  "\n"
)

alpha_true <- 1.1
stationary_var_true <- 0.9
ou_theta <- lapply(
  tree_list,
  simulate_ou_tips,
  alpha = alpha_true,
  stationary_var = stationary_var_true
)
ou_s <- lapply(ou_theta, function(theta) setNames(runif(length(theta), 0.1, 0.3), names(theta)))
ou_y_list <- Map(function(theta, s) {
  setNames(rnorm(length(theta), theta, s), names(theta))
}, ou_theta, ou_s)
ou_y <- unlist(ou_y_list)
ou_s_all <- unlist(ou_s)

ou_forest_fit <- ebnm_tree_ou(
  ou_y,
  ou_s_all,
  tree_list,
  g_init = tree_ou(alpha_true, stationary_var_true),
  fix_g = TRUE
)
ou_sep_fit <- fit_separately(tree_list, ou_y_list, ou_s, ebnm_tree_ou, tree_ou(alpha_true, stationary_var_true))

cat("\nOU forest check\n")
cat(
  "max abs diff posterior mean:",
  signif(max(abs(ou_forest_fit$posterior$mean - ou_sep_fit$posterior_mean)), 4),
  "\n"
)
cat(
  "max abs diff posterior sd:",
  signif(max(abs(ou_forest_fit$posterior$sd - ou_sep_fit$posterior_sd)), 4),
  "\n"
)
cat(
  "forest loglik - sum separate:",
  signif(as.numeric(ou_forest_fit$log_likelihood) - ou_sep_fit$loglik, 4),
  "\n"
)
