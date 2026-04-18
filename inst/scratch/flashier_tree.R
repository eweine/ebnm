library(ape)
library(flashier)

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", helpers = FALSE, quiet = TRUE)
} else {
  library(ebnm)
}

n_cells <- 1000
n_genes <- 250
noise_sd <- 1

set.seed(1)
tree <- ape::rtree(n_cells)
tree$edge.length <- pmax(tree$edge.length, 0.05)
ntip <- ape::Ntip(tree)
edge_pre <- ape::reorder.phylo(tree, "cladewise")$edge
elen_pre <- ape::reorder.phylo(tree, "cladewise")$edge.length
root_node <- setdiff(tree$edge[, 1], tree$edge[, 2])

simulate_tree_signal <- function(family = c("bm", "ou"),
                                 bm_var = 1,
                                 alpha = 1,
                                 stationary_var = 1) {
  family <- match.arg(family)
  x_node <- numeric(ntip + tree$Nnode)

  if (family == "bm") {
    x_node[root_node] <- 0
    for (i in seq_len(nrow(edge_pre))) {
      p <- edge_pre[i, 1]
      ch <- edge_pre[i, 2]
      x_node[ch] <- rnorm(1, x_node[p], sqrt(bm_var * elen_pre[i]))
    }
  } else {
    x_node[root_node] <- rnorm(1, 0, sqrt(stationary_var))
    for (i in seq_len(nrow(edge_pre))) {
      p <- edge_pre[i, 1]
      ch <- edge_pre[i, 2]
      a <- exp(-alpha * elen_pre[i])
      q <- stationary_var * (1 - a^2)
      x_node[ch] <- rnorm(1, a * x_node[p], sqrt(q))
    }
  }

  setNames(x_node[seq_len(ntip)], tree$tip.label)
}

simulate_factor_model <- function(l_true) {
  f_true <- rexp(n_genes)
  f_true[sample(n_genes, round(n_genes / 2))] <- 0
  mean_true <- tcrossprod(l_true, f_true)
  y <- mean_true + matrix(rnorm(n_cells * n_genes, sd = noise_sd), n_cells, n_genes)
  list(y = y, f_true = f_true)
}

fit_flash_tree <- function(y, tree_ebnm_fn) {
  flash(
    data = y,
    S = noise_sd,
    ebnm_fn = list(tree_ebnm_fn, ebnm::ebnm_point_exponential)
  )
}

align_factor <- function(est, truth) {
  sgn <- sign(sum(est * truth))
  if (sgn == 0) sgn <- 1
  est * sgn
}

report_fit <- function(label, flash_fit, l_true, f_true) {
  l_est <- align_factor(flash_fit$L_pm[, 1], l_true)
  f_est <- align_factor(flash_fit$F_pm[, 1], f_true)

  cat("\n", label, "\n", sep = "")
  cat("cor(L_est, L_true) =", round(cor(l_est, l_true), 3), "\n")
  cat("cor(F_est, F_true) =", round(cor(f_est, f_true), 3), "\n")

  invisible(list(l_est = l_est, f_est = f_est))
}

bm_var_true <- 1.0
l_true_bm <- simulate_tree_signal("bm", bm_var = bm_var_true)
sim_bm <- simulate_factor_model(l_true_bm)
flash_bm <- fit_flash_tree(sim_bm$y, ebnm::ebnm_tree_bm_fn(tree))
bm_diag <- report_fit("BM example", flash_bm, l_true_bm, sim_bm$f_true)
cat("estimated BM variance =", round(flash_bm$L_ghat[[1]]$bm_var, 3), "\n")
cat("true BM variance      =", bm_var_true, "\n")

par(mfrow = c(2, 2))
plot(
  l_true_bm,
  bm_diag$l_est,
  main = "BM: loadings",
  xlab = "True",
  ylab = "Estimated"
)
abline(0, 1, col = 2)
plot(
  sim_bm$f_true,
  bm_diag$f_est,
  main = "BM: factors",
  xlab = "True",
  ylab = "Estimated"
)
abline(0, 1, col = 2)

alpha_true <- 1.25
stationary_var_true <- 0.8
l_true_ou <- simulate_tree_signal(
  "ou",
  alpha = alpha_true,
  stationary_var = stationary_var_true
)
sim_ou <- simulate_factor_model(l_true_ou)
flash_ou <- fit_flash_tree(sim_ou$y, ebnm::ebnm_tree_ou_fn(tree))
ou_diag <- report_fit("OU example", flash_ou, l_true_ou, sim_ou$f_true)
cat("estimated OU alpha            =", round(flash_ou$L_ghat[[1]]$alpha, 3), "\n")
cat("true OU alpha                 =", alpha_true, "\n")
cat(
  "estimated OU stationary var   =",
  round(flash_ou$L_ghat[[1]]$stationary_var, 3),
  "\n"
)
cat("true OU stationary var        =", stationary_var_true, "\n")

plot(
  l_true_ou,
  ou_diag$l_est,
  main = "OU: loadings",
  xlab = "True",
  ylab = "Estimated"
)
abline(0, 1, col = 2)
plot(
  sim_ou$f_true,
  ou_diag$f_est,
  main = "OU: factors",
  xlab = "True",
  ylab = "Estimated"
)
abline(0, 1, col = 2)
