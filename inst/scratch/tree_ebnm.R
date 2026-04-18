# install.packages("ape")
library(ape)

# -----------------------------
# Helpers
# -----------------------------

.align_tree_data <- function(tree, y, s) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be an object of class 'phylo'.")
  }
  if (!is.rooted(tree)) {
    stop("tree must be rooted.")
  }
  if (!is.null(tree$root.edge)) {
    warning("tree$root.edge is ignored by these functions. If you want it included, absorb it into root_var manually.")
  }

  n <- length(tree$tip.label)

  if (length(y) != n || length(s) != n) {
    stop("y and s must each have length equal to the number of tips in tree.")
  }

  if (!is.null(names(y)) || !is.null(names(s))) {
    if (is.null(names(y)) || is.null(names(s))) {
      stop("If either y or s is named, both must be named.")
    }
    if (!all(tree$tip.label %in% names(y))) {
      stop("All tree tip labels must appear in names(y).")
    }
    if (!all(tree$tip.label %in% names(s))) {
      stop("All tree tip labels must appear in names(s).")
    }
    y <- y[tree$tip.label]
    s <- s[tree$tip.label]
  } else {
    names(y) <- tree$tip.label
    names(s) <- tree$tip.label
  }

  if (any(!is.finite(y))) stop("y contains non-finite values.")
  if (any(!is.finite(s))) stop("s contains non-finite values.")
  if (any(s <= 0)) stop("All s_i must be strictly positive.")

  list(
    tree = tree,
    y = as.numeric(y),
    s = as.numeric(s),
    tip.label = tree$tip.label
  )
}

.default_init_bm_var <- function(y, C) {
  denom <- mean(diag(C))
  denom <- max(denom, 1e-8)
  max(0.5 * stats::var(y) / denom, 1e-6)
}

.default_init_root_var <- function(y) {
  max(0.5 * stats::var(y), 1e-6)
}

.combine_gaussian_messages <- function(means, vars) {
  # Product_k N(x; means[k], vars[k]) = const * N(x; mean, var)
  if (length(means) != length(vars)) stop("means and vars must have same length.")
  if (length(means) == 0L) stop("Need at least one message to combine.")
  if (any(vars <= 0)) stop("All variances in .combine_gaussian_messages must be > 0.")

  if (length(means) == 1L) {
    return(list(mean = means[1], var = vars[1], log_const = 0))
  }

  m <- means[1]
  v <- vars[1]
  log_const <- 0

  for (k in 2:length(means)) {
    a <- means[k]
    b <- vars[k]

    # N(x;m,v) * N(x;a,b) = N(m;a,v+b) * N(x; m_new, v_new)
    log_const <- log_const + dnorm(m, mean = a, sd = sqrt(v + b), log = TRUE)
    v_new <- 1 / (1 / v + 1 / b)
    m_new <- v_new * (m / v + a / b)

    m <- m_new
    v <- v_new
  }

  list(mean = m, var = v, log_const = log_const)
}

.build_tree_index <- function(tree) {
  ntip <- Ntip(tree)
  nnode <- tree$Nnode
  total_nodes <- ntip + nnode

  edge <- tree$edge
  edge_len <- tree$edge.length

  root <- setdiff(edge[, 1], edge[, 2])
  if (length(root) != 1L) {
    stop("Could not identify a unique root.")
  }

  parent_of <- rep(NA_integer_, total_nodes)
  edge_len_to_child <- rep(NA_real_, total_nodes)
  children <- vector("list", total_nodes)

  for (i in seq_len(nrow(edge))) {
    p <- edge[i, 1]
    ch <- edge[i, 2]
    parent_of[ch] <- p
    edge_len_to_child[ch] <- edge_len[i]
    children[[p]] <- c(children[[p]], ch)
  }

  # Postorder internal-node order
  tree_post <- reorder.phylo(tree, order = "postorder")
  post_parents <- tree_post$edge[, 1]
  internal_post <- post_parents[!duplicated(post_parents, fromLast = TRUE)]

  # Preorder internal-node order
  tree_pre <- reorder.phylo(tree, order = "cladewise")
  internal_pre <- unique(tree_pre$edge[, 1])

  list(
    ntip = ntip,
    nnode = nnode,
    total_nodes = total_nodes,
    root = root,
    parent_of = parent_of,
    edge_len_to_child = edge_len_to_child,
    children = children,
    internal_post = internal_post,
    internal_pre = internal_pre
  )
}

# -----------------------------
# Dense solver
# -----------------------------

.prep_dense <- function(tree, y, s) {
  dat <- .align_tree_data(tree, y, s)

  C <- vcv.phylo(dat$tree)
  C <- C[dat$tip.label, dat$tip.label, drop = FALSE]

  list(
    tree = dat$tree,
    y = dat$y,
    s = dat$s,
    s2 = dat$s^2,
    C = C,
    n = length(dat$y),
    tip.label = dat$tip.label
  )
}

.eval_dense_given_par <- function(prep, bm_var, root_var) {
  if (bm_var < 0) stop("bm_var must be nonnegative.")
  if (root_var < 0) stop("root_var must be nonnegative.")

  n <- prep$n
  J <- matrix(1, n, n)
  K_theta <- bm_var * prep$C + root_var * J
  Sigma_y <- K_theta + diag(prep$s2, n, n)

  # Dense solve
  Sigma_y_inv_y <- solve(Sigma_y, prep$y)
  logdet <- as.numeric(determinant(Sigma_y, logarithm = TRUE)$modulus)

  loglik <- -0.5 * (
    n * log(2 * pi) +
      logdet +
      sum(prep$y * Sigma_y_inv_y)
  )

  post_mean <- drop(K_theta %*% Sigma_y_inv_y)
  post_cov <- K_theta - K_theta %*% solve(Sigma_y, K_theta)
  post_cov <- 0.5 * (post_cov + t(post_cov))

  names(post_mean) <- prep$tip.label
  names(diag(post_cov)) <- prep$tip.label

  list(
    loglik = loglik,
    bm_var = bm_var,
    root_var = root_var,
    prior_cov = K_theta,
    marginal_cov = Sigma_y,
    posterior_mean = post_mean,
    posterior_var = diag(post_cov),
    posterior_cov = post_cov
  )
}

eb_normalmeans_tree_dense <- function(
    tree,
    y,
    s,
    estimate_root_var = FALSE,
    root_var = 1,
    init_bm_var = NULL,
    init_root_var = NULL,
    lower_bm_var = 1e-10,
    lower_root_var = 1e-10,
    optim_control = list()
) {
  prep <- .prep_dense(tree, y, s)

  if (is.null(init_bm_var)) {
    init_bm_var <- .default_init_bm_var(prep$y, prep$C)
  }
  if (is.null(init_root_var)) {
    init_root_var <- .default_init_root_var(prep$y)
  }

  if (!estimate_root_var && root_var < 0) {
    stop("Fixed root_var must be nonnegative.")
  }

  if (estimate_root_var) {
    obj <- function(par_log) {
      bm_var_now <- exp(par_log[1])
      root_var_now <- exp(par_log[2])
      - .eval_dense_given_par(prep, bm_var_now, root_var_now)$loglik
    }

    opt <- optim(
      par = log(c(init_bm_var, init_root_var)),
      fn = obj,
      method = "L-BFGS-B",
      lower = log(c(lower_bm_var, lower_root_var)),
      control = optim_control
    )

    bm_var_hat <- exp(opt$par[1])
    root_var_hat <- exp(opt$par[2])
  } else {
    obj <- function(par_log) {
      bm_var_now <- exp(par_log[1])
      - .eval_dense_given_par(prep, bm_var_now, root_var)$loglik
    }

    opt <- optim(
      par = log(init_bm_var),
      fn = obj,
      method = "L-BFGS-B",
      lower = log(lower_bm_var),
      control = optim_control
    )

    bm_var_hat <- exp(opt$par[1])
    root_var_hat <- root_var
  }

  fit <- .eval_dense_given_par(prep, bm_var_hat, root_var_hat)

  fit$method <- "dense"
  fit$estimate_root_var <- estimate_root_var
  fit$optim <- opt
  fit
}

# -----------------------------
# Fast pruning / message-passing solver
# -----------------------------

.prep_pruning <- function(tree, y, s) {
  dat <- .align_tree_data(tree, y, s)
  idx <- .build_tree_index(dat$tree)

  list(
    tree = dat$tree,
    y = dat$y,
    s = dat$s,
    s2 = dat$s^2,
    tip.label = dat$tip.label,
    ntip = idx$ntip,
    nnode = idx$nnode,
    total_nodes = idx$total_nodes,
    root = idx$root,
    parent_of = idx$parent_of,
    edge_len_to_child = idx$edge_len_to_child,
    children = idx$children,
    internal_post = idx$internal_post,
    internal_pre = idx$internal_pre
  )
}

.upward_pass <- function(prep, bm_var) {
  if (bm_var < 0) stop("bm_var must be nonnegative.")

  total_nodes <- prep$total_nodes
  ntip <- prep$ntip

  up_mean <- rep(NA_real_, total_nodes)
  up_var <- rep(NA_real_, total_nodes)
  up_logconst <- rep(NA_real_, total_nodes)

  # Tip messages: p(y_i | x_i) as a function of x_i
  up_mean[1:ntip] <- prep$y
  up_var[1:ntip] <- prep$s2
  up_logconst[1:ntip] <- 0

  for (v in prep$internal_post) {
    ch <- prep$children[[v]]
    if (length(ch) == 0L) next

    msg_means <- up_mean[ch]
    msg_vars <- up_var[ch] + bm_var * prep$edge_len_to_child[ch]
    comb <- .combine_gaussian_messages(msg_means, msg_vars)

    up_mean[v] <- comb$mean
    up_var[v] <- comb$var
    up_logconst[v] <- sum(up_logconst[ch]) + comb$log_const
  }

  list(
    up_mean = up_mean,
    up_var = up_var,
    up_logconst = up_logconst
  )
}

.loglik_from_upward <- function(prep, up, root_var) {
  root <- prep$root
  m <- up$up_mean[root]
  v <- up$up_var[root]
  lc <- up$up_logconst[root]

  if (root_var < 0) stop("root_var must be nonnegative.")

  if (root_var == 0) {
    # root fixed at 0
    lc + dnorm(0, mean = m, sd = sqrt(v), log = TRUE)
  } else {
    lc + dnorm(m, mean = 0, sd = sqrt(v + root_var), log = TRUE)
  }
}

.downward_pass <- function(prep, up, bm_var, root_var) {
  total_nodes <- prep$total_nodes
  ntip <- prep$ntip
  root <- prep$root

  out_mean <- rep(NA_real_, total_nodes)
  out_var <- rep(NA_real_, total_nodes)

  # Outside message at the root is the root prior
  out_mean[root] <- 0
  out_var[root] <- root_var

  for (p in prep$internal_pre) {
    ch <- prep$children[[p]]
    if (length(ch) == 0L) next

    # Special case: root fixed exactly at 0
    if (p == root && root_var == 0) {
      for (c in ch) {
        out_mean[c] <- 0
        out_var[c] <- bm_var * prep$edge_len_to_child[c]
      }
      next
    }

    child_means_to_p <- up$up_mean[ch]
    child_vars_to_p <- up$up_var[ch] + bm_var * prep$edge_len_to_child[ch]

    for (c in ch) {
      sib <- ch[ch != c]

      excl_means <- c(out_mean[p], child_means_to_p[ch != c])
      excl_vars <- c(out_var[p], child_vars_to_p[ch != c])

      comb_excl <- .combine_gaussian_messages(excl_means, excl_vars)

      out_mean[c] <- comb_excl$mean
      out_var[c] <- comb_excl$var + bm_var * prep$edge_len_to_child[c]
    }
  }

  # Posterior marginals for all nodes
  post_mean_node <- rep(NA_real_, total_nodes)
  post_var_node <- rep(NA_real_, total_nodes)

  # Root
  if (root_var == 0) {
    post_mean_node[root] <- 0
    post_var_node[root] <- 0
  } else {
    vr <- 1 / (1 / root_var + 1 / up$up_var[root])
    mr <- vr * (0 / root_var + up$up_mean[root] / up$up_var[root])
    post_mean_node[root] <- mr
    post_var_node[root] <- vr
  }

  # Non-root nodes
  for (v in setdiff(seq_len(total_nodes), root)) {
    pv <- 1 / (1 / out_var[v] + 1 / up$up_var[v])
    pm <- pv * (out_mean[v] / out_var[v] + up$up_mean[v] / up$up_var[v])
    post_mean_node[v] <- pm
    post_var_node[v] <- pv
  }

  list(
    out_mean = out_mean,
    out_var = out_var,
    posterior_node_mean = post_mean_node,
    posterior_node_var = post_var_node,
    posterior_tip_mean = stats::setNames(post_mean_node[1:ntip], prep$tip.label),
    posterior_tip_var = stats::setNames(post_var_node[1:ntip], prep$tip.label)
  )
}

.eval_pruning_given_par <- function(prep, bm_var, root_var) {
  up <- .upward_pass(prep, bm_var)
  loglik <- .loglik_from_upward(prep, up, root_var)
  down <- .downward_pass(prep, up, bm_var, root_var)

  list(
    loglik = loglik,
    bm_var = bm_var,
    root_var = root_var,
    posterior_mean = down$posterior_tip_mean,
    posterior_var = down$posterior_tip_var,
    posterior_node_mean = down$posterior_node_mean,
    posterior_node_var = down$posterior_node_var,
    upward = up,
    downward = down
  )
}

eb_normalmeans_tree_pruning <- function(
    tree,
    y,
    s,
    estimate_root_var = FALSE,
    root_var = 1,
    init_bm_var = NULL,
    init_root_var = NULL,
    lower_bm_var = 1e-10,
    lower_root_var = 1e-10,
    optim_control = list()
) {
  prep <- .prep_pruning(tree, y, s)

  # For initial values, use the dense covariance diagonal via vcv.phylo just for a rough scale.
  C_for_init <- vcv.phylo(prep$tree)
  C_for_init <- C_for_init[prep$tip.label, prep$tip.label, drop = FALSE]

  if (is.null(init_bm_var)) {
    init_bm_var <- .default_init_bm_var(prep$y, C_for_init)
  }
  if (is.null(init_root_var)) {
    init_root_var <- .default_init_root_var(prep$y)
  }

  if (!estimate_root_var && root_var < 0) {
    stop("Fixed root_var must be nonnegative.")
  }

  if (estimate_root_var) {
    obj <- function(par_log) {
      bm_var_now <- exp(par_log[1])
      root_var_now <- exp(par_log[2])
      - .eval_pruning_given_par(prep, bm_var_now, root_var_now)$loglik
    }

    opt <- optim(
      par = log(c(init_bm_var, init_root_var)),
      fn = obj,
      method = "L-BFGS-B",
      lower = log(c(lower_bm_var, lower_root_var)),
      control = optim_control
    )

    bm_var_hat <- exp(opt$par[1])
    root_var_hat <- exp(opt$par[2])
  } else {
    obj <- function(par_log) {
      bm_var_now <- exp(par_log[1])
      - .eval_pruning_given_par(prep, bm_var_now, root_var)$loglik
    }

    opt <- optim(
      par = log(init_bm_var),
      fn = obj,
      method = "L-BFGS-B",
      lower = log(lower_bm_var),
      control = optim_control
    )

    bm_var_hat <- exp(opt$par[1])
    root_var_hat <- root_var
  }

  fit <- .eval_pruning_given_par(prep, bm_var_hat, root_var_hat)

  fit$method <- "pruning"
  fit$estimate_root_var <- estimate_root_var
  fit$optim <- opt
  fit
}

set.seed(123)

# -----------------------------
# 1. Simulate a random tree
# -----------------------------
n_tips <- 800
tree <- rtree(n_tips)

# Make sure branch lengths are all positive and not too tiny
tree$edge.length <- pmax(tree$edge.length, 0.05)

plot(tree, main = "Simulated tree")

# -----------------------------
# 2. Simulate latent states on the tree
# Model:
#   x_root ~ N(0, root_var_true)
#   x_child | x_parent ~ N(x_parent, bm_var_true * branch_length)
# -----------------------------
bm_var_true   <- 0.1
root_var_true <- 3

simulate_bm_tree_states <- function(tree, bm_var, root_var, root_mean = 0) {
  ntip  <- Ntip(tree)
  nnode <- tree$Nnode
  total_nodes <- ntip + nnode

  edge <- tree$edge
  elen <- tree$edge.length

  root <- setdiff(edge[, 1], edge[, 2])
  if (length(root) != 1L) stop("Tree must have a unique root.")

  # reorder so parents appear before children
  tree_pre <- reorder.phylo(tree, order = "cladewise")
  edge <- tree_pre$edge
  elen <- tree_pre$edge.length

  x <- rep(NA_real_, total_nodes)
  x[root] <- rnorm(1, mean = root_mean, sd = sqrt(root_var))

  for (k in seq_len(nrow(edge))) {
    p <- edge[k, 1]
    c <- edge[k, 2]
    tpc <- elen[k]
    x[c] <- rnorm(1, mean = x[p], sd = sqrt(bm_var * tpc))
  }

  list(
    node_states = x,
    tip_states = stats::setNames(x[1:ntip], tree$tip.label),
    root = root
  )
}

sim_states <- simulate_bm_tree_states(
  tree = tree,
  bm_var = bm_var_true,
  root_var = root_var_true
)

theta_true <- sim_states$tip_states

# -----------------------------
# 3. Add observation noise
# y_i ~ N(theta_i, s_i^2)
# -----------------------------
s <- runif(n_tips, min = 0.15, max = 0.5)
names(s) <- tree$tip.label

y <- rnorm(n_tips, mean = theta_true, sd = s)
names(y) <- tree$tip.label

# Quick look
plot(theta_true, y,
     xlab = expression(true~theta[i]),
     ylab = expression(observed~y[i]),
     main = "Observed values vs true latent tip states")
abline(0, 1, col = "red", lty = 2)

# -----------------------------
# 4. Fit EB solvers
#    (a) fixed root variance
# -----------------------------
fit_dense_fixed <- eb_normalmeans_tree_dense(
  tree = tree,
  y = y,
  s = s,
  estimate_root_var = FALSE,
  root_var = 1000
)

fit_fast_fixed <- eb_normalmeans_tree_pruning(
  tree = tree,
  y = y,
  s = s,
  estimate_root_var = FALSE,
  root_var = 1000
)

# -----------------------------
# 5. Fit EB solvers
#    (b) estimate root variance too
# -----------------------------
fit_dense_est <- eb_normalmeans_tree_dense(
  tree = tree,
  y = y,
  s = s,
  estimate_root_var = TRUE
)

fit_fast_est <- eb_normalmeans_tree_pruning(
  tree = tree,
  y = y,
  s = s,
  estimate_root_var = TRUE
)

# -----------------------------
# 6. Compare results
# -----------------------------
cat("\nTrue hyperparameters:\n")
cat("  bm_var   =", bm_var_true, "\n")
cat("  root_var =", root_var_true, "\n")

cat("\nDense, fixed root variance:\n")
cat("  bm_var   =", fit_dense_fixed$bm_var, "\n")
cat("  root_var =", fit_dense_fixed$root_var, "\n")
cat("  loglik   =", fit_dense_fixed$loglik, "\n")

cat("\nPruning, fixed root variance:\n")
cat("  bm_var   =", fit_fast_fixed$bm_var, "\n")
cat("  root_var =", fit_fast_fixed$root_var, "\n")
cat("  loglik   =", fit_fast_fixed$loglik, "\n")

cat("\nDense, estimated root variance:\n")
cat("  bm_var   =", fit_dense_est$bm_var, "\n")
cat("  root_var =", fit_dense_est$root_var, "\n")
cat("  loglik   =", fit_dense_est$loglik, "\n")

cat("\nPruning, estimated root variance:\n")
cat("  bm_var   =", fit_fast_est$bm_var, "\n")
cat("  root_var =", fit_fast_est$root_var, "\n")
cat("  loglik   =", fit_fast_est$loglik, "\n")

# Dense and pruning should agree closely
cat("\nMax abs difference in posterior means (estimated-root fits):\n")
cat(max(abs(fit_dense_est$posterior_mean - fit_fast_est$posterior_mean)), "\n")

cat("\nMax abs difference in posterior variances (estimated-root fits):\n")
cat(max(abs(fit_dense_est$posterior_var - fit_fast_est$posterior_var)), "\n")

# -----------------------------
# 7. Compare posterior means to truth
# -----------------------------
op <- par(mfrow = c(1, 2))

plot(theta_true, fit_dense_est$posterior_mean,
     xlab = expression(true~theta[i]),
     ylab = "Posterior mean",
     main = "Dense EB posterior mean")
abline(0, 1, col = "red", lty = 2)

plot(theta_true, fit_fast_est$posterior_mean,
     xlab = expression(true~theta[i]),
     ylab = "Posterior mean",
     main = "Pruning EB posterior mean")
abline(0, 1, col = "red", lty = 2)

par(op)

# -----------------------------
# 8. Simple summary metrics
# -----------------------------
rmse_obs <- sqrt(mean((y - theta_true)^2))
rmse_post_dense <- sqrt(mean((fit_dense_est$posterior_mean - theta_true)^2))
rmse_post_fast  <- sqrt(mean((fit_fast_est$posterior_mean - theta_true)^2))

cat("\nRMSE of raw observations vs truth: ", rmse_obs, "\n")
cat("RMSE of dense posterior mean vs truth: ", rmse_post_dense, "\n")
cat("RMSE of pruning posterior mean vs truth: ", rmse_post_fast, "\n")
