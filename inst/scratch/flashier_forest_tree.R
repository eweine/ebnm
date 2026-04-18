library(ape)
library(flashier)

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", helpers = FALSE, quiet = TRUE)
} else {
  library(ebnm)
}

set.seed(1)

n_genes <- 250
noise_sd <- 1

make_tree <- function(n_tip, prefix) {
  tr <- ape::rtree(n_tip, tip.label = paste0(prefix, seq_len(n_tip)))
  tr$edge.length <- pmax(tr$edge.length, 0.05)
  tr
}

tree_list <- list(
  make_tree(300, "a_"),
  make_tree(350, "b_"),
  make_tree(250, "c_")
)
n_cells <- sum(vapply(tree_list, ape::Ntip, integer(1)))

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

align_factor <- function(est, truth) {
  sgn <- sign(sum(est * truth))
  if (sgn == 0) sgn <- 1
  est * sgn
}

simulate_factor_model <- function(l_true) {
  f_true <- rexp(n_genes)
  f_true[sample(n_genes, round(n_genes / 2))] <- 0
  mean_true <- tcrossprod(l_true, f_true)
  y <- mean_true + matrix(rnorm(n_cells * n_genes, sd = noise_sd), n_cells, n_genes)
  list(y = y, f_true = f_true)
}

fit_flash_forest <- function(y, tree_ebnm_fn) {
  flash(
    data = y,
    S = noise_sd,
    ebnm_fn = list(tree_ebnm_fn, ebnm::ebnm_point_exponential)
  )
}

report_fit <- function(label, flash_fit, l_true, f_true) {
  l_est <- align_factor(flash_fit$L_pm[, 1], l_true)
  f_est <- align_factor(flash_fit$F_pm[, 1], f_true)
  slope_l <- unname(stats::coef(stats::lm(l_est ~ 0 + l_true)))

  cat("\n", label, "\n", sep = "")
  cat("forest sizes              =", paste(vapply(tree_list, ape::Ntip, integer(1)), collapse = ", "), "\n")
  cat("cor(L_est, L_true)        =", round(stats::cor(l_est, l_true), 3), "\n")
  cat("cor(F_est, F_true)        =", round(stats::cor(f_est, f_true), 3), "\n")
  cat("loading rescaling slope   =", round(slope_l, 3), "\n")

  invisible(list(l_est = l_est, f_est = f_est, slope_l = slope_l))
}

bm_var_true <- 1
l_true_bm <- unlist(lapply(tree_list, simulate_bm_tips, bm_var = bm_var_true))
sim_bm <- simulate_factor_model(l_true_bm)
flash_bm <- fit_flash_forest(sim_bm$y, ebnm::ebnm_tree_bm_fn(tree_list))
bm_diag <- report_fit("BM forest example", flash_bm, l_true_bm, sim_bm$f_true)
cat("estimated BM variance     =", round(flash_bm$L_ghat[[1]]$bm_var, 3), "\n")
cat("true BM variance          =", bm_var_true, "\n")
cat("scaled true BM variance   =", round(bm_var_true * bm_diag$slope_l^2, 3), "\n")

par(mfrow = c(2, 2))
plot(
  l_true_bm,
  bm_diag$l_est,
  main = "BM forest: loadings",
  xlab = "True",
  ylab = "Estimated"
)
abline(0, 1, col = 2)
plot(
  sim_bm$f_true,
  bm_diag$f_est,
  main = "BM forest: factors",
  xlab = "True",
  ylab = "Estimated"
)
abline(0, 1, col = 2)

alpha_true <- 1.2
stationary_var_true <- 0.8
l_true_ou <- unlist(lapply(
  tree_list,
  simulate_ou_tips,
  alpha = alpha_true,
  stationary_var = stationary_var_true
))
sim_ou <- simulate_factor_model(l_true_ou)
flash_ou <- fit_flash_forest(sim_ou$y, ebnm::ebnm_tree_ou_fn(tree_list))
ou_diag <- report_fit("OU forest example", flash_ou, l_true_ou, sim_ou$f_true)
cat("estimated OU alpha        =", round(flash_ou$L_ghat[[1]]$alpha, 3), "\n")
cat("true OU alpha             =", alpha_true, "\n")
cat("estimated OU stat var     =", round(flash_ou$L_ghat[[1]]$stationary_var, 3), "\n")
cat("true OU stat var          =", stationary_var_true, "\n")
cat(
  "scaled true OU stat var   =",
  round(stationary_var_true * ou_diag$slope_l^2, 3),
  "\n"
)

plot(
  l_true_ou,
  ou_diag$l_est,
  main = "OU forest: loadings",
  xlab = "True",
  ylab = "Estimated"
)
abline(0, 1, col = 2)
plot(
  sim_ou$f_true,
  ou_diag$f_est,
  main = "OU forest: factors",
  xlab = "True",
  ylab = "Estimated"
)
abline(0, 1, col = 2)
