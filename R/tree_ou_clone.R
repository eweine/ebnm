#' Constructor for tree_ou_clone prior class
#'
#' Creates a hierarchical OU forest prior with shared OU parameters and
#' tree-specific means drawn from a mean-zero normal distribution.
#'
#' @param alpha A positive scalar giving the OU reversion rate per unit
#'   branch length.
#' @param stationary_var A positive scalar giving the stationary variance of
#'   the OU process around each tree-specific mean.
#' @param mean_var A positive scalar giving the variance of the tree-specific
#'   means.
#'
#' @return An object of class \code{tree_ou_clone}.
#'
#' @export
tree_ou_clone <- function(alpha, stationary_var, mean_var) {
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0) {
    stop("alpha must be a single positive number.")
  }
  if (!is.numeric(stationary_var) || length(stationary_var) != 1L || stationary_var <= 0) {
    stop("stationary_var must be a single positive number.")
  }
  if (!is.numeric(mean_var) || length(mean_var) != 1L || mean_var <= 0) {
    stop("mean_var must be a single positive number.")
  }
  structure(
    list(alpha = alpha, stationary_var = stationary_var, mean_var = mean_var),
    class = "tree_ou_clone"
  )
}

#' Create a flashier-compatible ebnm function for the tree OU clone prior
#'
#' @param tree A rooted phylogenetic tree of class \code{\link[ape]{phylo}},
#'   or a list of rooted \code{phylo} objects representing independent trees.
#' @param tree_index An optional pre-computed index from
#'   \code{\link{tree_ou_clone_precomp}}. If \code{NULL}, it is built
#'   automatically from \code{tree}.
#'
#' @return A function suitable for use as the \code{ebnm_fn} argument to
#'   \code{\link[flashier]{flash}} or \code{\link[flashier]{flash_greedy}}.
#'
#' @export
ebnm_tree_ou_clone_fn <- function(tree, tree_index = NULL) {
  if (is.null(tree_index)) {
    tree_index <- tree_ou_clone_precomp(tree)
  }
  ntip <- tree_index$ntip
  function(x, s, g_init = NULL, fix_g = FALSE, output) {
    if (length(x) != ntip) {
      return(ebnm_normal(x, s, g_init = NULL, fix_g = FALSE, output = output))
    }
    ebnm_tree_ou_clone(
      x = x,
      s = s,
      tree = tree,
      tree_index = tree_index,
      g_init = g_init,
      fix_g = fix_g,
      output = output
    )
  }
}

#' Pre-compute tree structure for repeated ebnm_tree_ou_clone calls
#'
#' @param tree A rooted phylogenetic tree of class \code{\link[ape]{phylo}},
#'   or a list of rooted \code{phylo} objects representing independent trees.
#'
#' @return A list of class \code{tree_ou_clone_index}.
#'
#' @export
tree_ou_clone_precomp <- function(tree) {
  idx <- .tree_phylo_precomp(tree)
  structure(idx, class = c("tree_ou_clone_index", class(idx)))
}

#' Solve the EBNM problem using a hierarchical OU forest prior
#'
#' Fits an empirical Bayes normal means model on a rooted tree or forest in
#' which each tree has its own latent mean \eqn{\mu_t}, with
#' \eqn{\mu_t \sim N(0, \sigma_\mu^2)}. Conditional on \eqn{\mu_t}, each tree
#' follows a stationary OU process centered at \eqn{\mu_t} with shared
#' reversion rate \eqn{\alpha} and stationary variance \eqn{\tau^2}.
#'
#' @param x A numeric vector of observations at tree tips. If named, names
#'   must match the concatenated tip labels from \code{tree}; output is always
#'   in that order.
#' @param s A scalar standard error applied to all tips, or a vector of
#'   per-tip standard errors of the same length as \code{x}. Must be strictly
#'   positive.
#' @param tree A rooted phylogenetic tree of class \code{\link[ape]{phylo}},
#'   or a list of rooted \code{phylo} objects. Different trees are modeled as
#'   prior-independent given their tree-specific means.
#' @param tree_index An optional pre-computed index returned by
#'   \code{\link{tree_ou_clone_precomp}}.
#' @param g_init An optional \code{tree_ou_clone} object (or an \code{ebnm}
#'   object whose \code{fitted_g} is a \code{tree_ou_clone}) specifying an
#'   initial or fixed prior.
#' @param fix_g If \code{TRUE}, fix the prior at \code{g_init} without further
#'   optimization.
#' @param output A character vector of values to return.
#' @param control A named list of control parameters passed to
#'   \code{\link[stats]{optim}} during optimization.
#'
#' @return An \code{ebnm} object whose \code{fitted_g} element is a
#'   \code{tree_ou_clone} object with fields \code{alpha},
#'   \code{stationary_var}, and \code{mean_var}.
#'
#' @export
ebnm_tree_ou_clone <- function(
    x,
    s = 1,
    tree,
    tree_index = NULL,
    g_init = NULL,
    fix_g = FALSE,
    output = ebnm_output_default(),
    control = NULL
) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' must be installed to use ebnm_tree_ou_clone.")
  }

  call <- match.call()
  dat <- .tbm_align_data(tree, x, s)
  y <- dat$y
  sv <- dat$s
  s2 <- sv^2
  tip.label <- dat$tip.label

  idx <- .tou_clone_validate_tree_index(tree, tree_index)

  if (inherits(g_init, "ebnm")) {
    g_init <- g_init[["fitted_g"]]
  }
  if (!is.null(g_init) && !inherits(g_init, "tree_ou_clone")) {
    stop("g_init must be an object of class 'tree_ou_clone'.")
  }
  if (fix_g && is.null(g_init)) {
    stop("g_init must be provided when fix_g = TRUE.")
  }

  if (fix_g) {
    alpha_hat <- g_init$alpha
    stationary_var_hat <- g_init$stationary_var
    mean_var_hat <- g_init$mean_var
  } else {
    bounds <- .tou_clone_bounds(y, sv, idx, control)
    par_init <- .tou_clone_init_par(y, sv, idx, g_init, bounds$lower, bounds$upper)

    negloglik_fn <- function(par) {
      -.tou_clone_fit_forest(
        y = y,
        s2 = s2,
        idx = idx,
        alpha = par[[1]],
        stationary_var = par[[2]],
        mean_var = par[[3]],
        compute_posterior = FALSE
      )$loglik
    }

    optim_control <- control
    optim_control$par_init <- NULL
    optim_control$lower <- NULL
    optim_control$upper <- NULL

    opt <- stats::optim(
      par = par_init,
      fn = negloglik_fn,
      method = "L-BFGS-B",
      lower = bounds$lower,
      upper = bounds$upper,
      control = optim_control
    )

    alpha_hat <- opt$par[[1]]
    stationary_var_hat <- opt$par[[2]]
    mean_var_hat <- opt$par[[3]]
  }

  fit <- .tou_clone_fit_forest(
    y = y,
    s2 = s2,
    idx = idx,
    alpha = alpha_hat,
    stationary_var = stationary_var_hat,
    mean_var = mean_var_hat,
    compute_posterior = TRUE
  )

  post_mean <- stats::setNames(fit$post_tip_mean, tip.label)
  post_var <- stats::setNames(fit$post_tip_var, tip.label)
  post_sd <- sqrt(post_var)
  y_named <- stats::setNames(y, tip.label)
  s_named <- stats::setNames(sv, tip.label)

  retlist <- list()
  if (data_in_output(output)) {
    retlist <- add_data_to_retlist(retlist, y_named, s_named)
  }
  if (posterior_in_output(output)) {
    posterior <- list()
    if (result_in_output(output)) {
      posterior$mean <- post_mean
      posterior$sd <- post_sd
      posterior$mean2 <- post_mean^2 + post_var
    }
    if (lfsr_in_output(output)) {
      posterior$lfsr <- stats::pnorm(0, abs(post_mean), post_sd)
    }
    retlist <- add_posterior_to_retlist(retlist, posterior, output, y_named)
  }
  if (g_in_output(output)) {
    retlist <- add_g_to_retlist(
      retlist,
      tree_ou_clone(alpha_hat, stationary_var_hat, mean_var_hat)
    )
  }
  if (llik_in_output(output)) {
    df <- if (fix_g) 0L else 3L
    retlist <- add_llik_to_retlist(retlist, fit$loglik, y_named, df = df)
  }
  if (sampler_in_output(output)) {
    warning("Posterior sampler is not implemented for ebnm_tree_ou_clone.")
  }

  as_ebnm(retlist, call)
}

.tou_clone_validate_tree_index <- function(tree, tree_index) {
  if (is.null(tree_index)) {
    return(tree_ou_clone_precomp(tree))
  }
  if (!inherits(tree_index, "tree_ou_clone_index") &&
      !inherits(tree_index, "tree_ou_index") &&
      !inherits(tree_index, "tree_bm_index")) {
    stop("tree_index must be an object returned by tree_ou_clone_precomp().")
  }
  tree_index
}

.tou_clone_bounds <- function(y, s, idx, control) {
  min_edge <- max(min(idx$edge_len_to_child[!is.na(idx$edge_len_to_child)]), 1e-8)
  scale_upper <- max(stats::var(y), stats::median(s^2), 1.0) * 100
  default_lower <- c(alpha = 1e-8, stationary_var = 1e-10, mean_var = 1e-10)
  default_upper <- c(alpha = 100 / min_edge, stationary_var = scale_upper, mean_var = scale_upper)

  lower <- default_lower
  upper <- default_upper
  if (!is.null(control$lower)) {
    lower <- rep_len(control$lower, 3)
    names(lower) <- names(default_lower)
  }
  if (!is.null(control$upper)) {
    upper <- rep_len(control$upper, 3)
    names(upper) <- names(default_upper)
  }
  list(lower = lower, upper = upper)
}

.tou_clone_init_par <- function(y, s, idx, g_init, lower, upper) {
  if (!is.null(g_init)) {
    par <- c(
      alpha = g_init$alpha,
      stationary_var = g_init$stationary_var,
      mean_var = g_init$mean_var
    )
  } else {
    edge_lengths <- idx$edge_len_to_child[!is.na(idx$edge_len_to_child)]
    signal_var <- max(stats::var(y) - stats::median(s^2), 1e-4)
    par <- c(
      alpha = 1 / max(stats::median(edge_lengths), 1e-6),
      stationary_var = signal_var / 2,
      mean_var = signal_var / 2
    )
  }
  pmin(pmax(par, lower), upper)
}

.tou_clone_fit_forest <- function(y, s2, idx, alpha, stationary_var, mean_var, compute_posterior) {
  delta <- 1
  post_mean <- numeric(idx$ntip)
  post_var <- numeric(idx$ntip)
  loglik <- 0

  for (k in seq_len(idx$n_tree)) {
    tip_idx <- idx$tree_tip_idx[[k]]
    subidx <- idx$subindex[[k]]
    yk <- y[tip_idx]
    s2k <- s2[tip_idx]

    fit0 <- tree_ou_compute_cpp(
      ntip = subidx$ntip,
      total_nodes = subidx$total_nodes,
      roots = subidx$root,
      y = yk,
      s2 = s2k,
      edge_len_to_child = subidx$edge_len_to_child,
      children = subidx$children,
      internal_post = subidx$internal_post,
      internal_pre = subidx$internal_pre,
      alpha = alpha,
      stationary_var = stationary_var
    )
    fit_pos <- tree_ou_compute_cpp(
      ntip = subidx$ntip,
      total_nodes = subidx$total_nodes,
      roots = subidx$root,
      y = yk - delta,
      s2 = s2k,
      edge_len_to_child = subidx$edge_len_to_child,
      children = subidx$children,
      internal_post = subidx$internal_post,
      internal_pre = subidx$internal_pre,
      alpha = alpha,
      stationary_var = stationary_var
    )
    fit_neg <- tree_ou_compute_cpp(
      ntip = subidx$ntip,
      total_nodes = subidx$total_nodes,
      roots = subidx$root,
      y = yk + delta,
      s2 = s2k,
      edge_len_to_child = subidx$edge_len_to_child,
      children = subidx$children,
      internal_post = subidx$internal_post,
      internal_pre = subidx$internal_pre,
      alpha = alpha,
      stationary_var = stationary_var
    )

    quad_a <- (fit_pos$loglik + fit_neg$loglik - 2 * fit0$loglik) / (2 * delta^2)
    quad_b <- (fit_pos$loglik - fit_neg$loglik) / (2 * delta)
    quad_c <- fit0$loglik

    precision_post <- 1 / mean_var - 2 * quad_a
    mean_post <- quad_b / precision_post
    var_post <- 1 / precision_post

    loglik <- loglik +
      quad_c -
      0.5 * log(mean_var) -
      0.5 * log(precision_post) +
      0.5 * quad_b^2 / precision_post

    if (compute_posterior) {
      cond0 <- fit0$post_tip_mean
      cond1 <- delta + fit_pos$post_tip_mean
      slope <- (cond1 - cond0) / delta
      post_mean[tip_idx] <- cond0 + slope * mean_post
      post_var[tip_idx] <- fit0$post_tip_var + slope^2 * var_post
    }
  }

  list(loglik = loglik, post_tip_mean = post_mean, post_tip_var = post_var)
}
