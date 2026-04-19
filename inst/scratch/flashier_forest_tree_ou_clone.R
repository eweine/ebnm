library(ape)
library(flashier)

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", helpers = FALSE, quiet = TRUE)
} else {
  library(ebnm)
}

set.seed(1)

n_genes <- 300
noise_sd <- 1

make_tree <- function(n_tip, prefix) {
  tr <- ape::rtree(n_tip, tip.label = paste0(prefix, seq_len(n_tip)))
  tr$edge.length <- pmax(tr$edge.length, 0.05)
  tr
}

tree_list <- list(
  make_tree(350, "a_"),
  make_tree(325, "b_"),
  make_tree(375, "c_"),
  make_tree(300, "d_")
)
n_cells <- sum(vapply(tree_list, ape::Ntip, integer(1)))

simulate_ou_clone_tips <- function(tree, alpha, stationary_var, mean_var) {
  ntip <- ape::Ntip(tree)
  edge_pre <- ape::reorder.phylo(tree, "cladewise")$edge
  elen_pre <- ape::reorder.phylo(tree, "cladewise")$edge.length
  root_node <- setdiff(tree$edge[, 1], tree$edge[, 2])
  mu_tree <- rnorm(1, 0, sqrt(mean_var))
  x_node <- numeric(ntip + tree$Nnode)
  x_node[root_node] <- rnorm(1, mu_tree, sqrt(stationary_var))

  for (i in seq_len(nrow(edge_pre))) {
    p <- edge_pre[i, 1]
    ch <- edge_pre[i, 2]
    a <- exp(-alpha * elen_pre[i])
    q <- stationary_var * (1 - a^2)
    x_node[ch] <- rnorm(1, mu_tree + a * (x_node[p] - mu_tree), sqrt(q))
  }

  list(
    mu_tree = mu_tree,
    theta = setNames(x_node[seq_len(ntip)], tree$tip.label)
  )
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

alpha_true <- 1.1
stationary_var_true <- 0.75
mean_var_true <- 0.6

sim_list <- lapply(
  tree_list,
  simulate_ou_clone_tips,
  alpha = alpha_true,
  stationary_var = stationary_var_true,
  mean_var = mean_var_true
)

l_true <- unlist(lapply(sim_list, `[[`, "theta"))
mu_true <- vapply(sim_list, `[[`, numeric(1), "mu_tree")
sim <- simulate_factor_model(l_true)

flash_fit <- fit_flash_forest(sim$y, ebnm::ebnm_tree_ou_clone_fn(tree_list))

l_est <- align_factor(flash_fit$L_pm[, 1], l_true)
f_est <- align_factor(flash_fit$F_pm[, 1], sim$f_true)
slope_l <- unname(stats::coef(stats::lm(l_est ~ 0 + l_true)))
scaled_stationary_var <- stationary_var_true * slope_l^2
scaled_mean_var <- mean_var_true * slope_l^2

cat("OU clone forest example\n")
cat("forest sizes                    =", paste(vapply(tree_list, ape::Ntip, integer(1)), collapse = ", "), "\n")
cat("number of cells                 =", n_cells, "\n")
cat("number of genes                 =", n_genes, "\n")
cat("cor(L_est, L_true)              =", round(stats::cor(l_est, l_true), 3), "\n")
cat("cor(F_est, F_true)              =", round(stats::cor(f_est, sim$f_true), 3), "\n")
cat("loading rescaling slope         =", round(slope_l, 3), "\n")
cat("estimated alpha                 =", round(flash_fit$L_ghat[[1]]$alpha, 3), "\n")
cat("true alpha                      =", alpha_true, "\n")
cat("estimated stationary var        =", round(flash_fit$L_ghat[[1]]$stationary_var, 3), "\n")
cat("true stationary var             =", stationary_var_true, "\n")
cat("scaled true stationary var      =", round(scaled_stationary_var, 3), "\n")
cat("estimated tree-mean var         =", round(flash_fit$L_ghat[[1]]$mean_var, 3), "\n")
cat("true tree-mean var              =", mean_var_true, "\n")
cat("scaled true tree-mean var       =", round(scaled_mean_var, 3), "\n")
cat("realized tree means             =", paste(round(mu_true, 3), collapse = ", "), "\n")

par(mfrow = c(1, 2))
plot(
  l_true,
  l_est,
  main = "OU Clone Forest: loadings",
  xlab = "True",
  ylab = "Estimated"
)
abline(0, 1, col = 2)

plot(
  sim$f_true,
  f_est,
  main = "OU Clone Forest: factors",
  xlab = "True",
  ylab = "Estimated"
)
abline(0, 1, col = 2)
