#' Constructor for tree_ou prior class
#'
#' Creates an Ornstein-Uhlenbeck tree prior object storing the estimated
#' reversion rate and stationary variance.
#'
#' @param alpha A positive scalar giving the OU reversion rate per unit
#'   branch length.
#' @param stationary_var A positive scalar giving the stationary variance of
#'   the OU process.
#'
#' @return An object of class \code{tree_ou}.
#'
#' @export
tree_ou <- function(alpha, stationary_var) {
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("alpha must be a single positive number.")
  }
  if (!is.numeric(stationary_var) ||
      length(stationary_var) != 1 ||
      stationary_var <= 0) {
    stop("stationary_var must be a single positive number.")
  }
  structure(
    list(alpha = alpha, stationary_var = stationary_var),
    class = "tree_ou"
  )
}

#' Create a flashier-compatible ebnm function for the tree OU prior
#'
#' Returns a closure with the signature \code{function(x, s, g_init, fix_g,
#' output)} expected by \pkg{flashier}. The \code{tree} and the pre-computed
#' \code{tree_index} are captured in the closure, so they are built only once
#' even when the function is called many times inside an iterative algorithm.
#'
#' @param tree A rooted phylogenetic tree of class \code{\link[ape]{phylo}}.
#' @param tree_index An optional pre-computed index from
#'   \code{\link{tree_ou_precomp}}. If \code{NULL} (default), it is built
#'   automatically from \code{tree}.
#'
#' @return A function suitable for use as the \code{ebnm_fn} argument to
#'   \code{\link[flashier]{flash}} or \code{\link[flashier]{flash_greedy}}.
#'
#' @seealso \code{\link{ebnm_tree_ou}}, \code{\link{tree_ou_precomp}}
#'
#' @export
ebnm_tree_ou_fn <- function(tree, tree_index = NULL) {
  if (is.null(tree_index)) {
    tree_index <- tree_ou_precomp(tree)
  }
  ntip <- tree_index$ntip
  function(x, s, g_init = NULL, fix_g = FALSE, output) {
    if (length(x) != ntip) {
      return(ebnm_normal(x, s, g_init = NULL, fix_g = FALSE, output = output))
    }
    ebnm_tree_ou(
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

#' Pre-compute tree structure for repeated ebnm_tree_ou calls
#'
#' Builds the message-passing index for a phylogenetic tree. Passing this
#' pre-computed index to \code{\link{ebnm_tree_ou}} via the \code{tree_index}
#' argument avoids rebuilding it on every call, which is important when
#' \code{ebnm_tree_ou} is used iteratively (e.g., as a prior in an
#' alternating optimization).
#'
#' @param tree A rooted phylogenetic tree of class \code{\link[ape]{phylo}}.
#'
#' @return A list of class \code{tree_ou_index} containing the pre-computed
#'   tree structure needed by \code{\link{ebnm_tree_ou}}.
#'
#' @seealso \code{\link{ebnm_tree_ou}}
#'
#' @export
tree_ou_precomp <- function(tree) {
  idx <- .tree_phylo_precomp(tree)
  structure(idx, class = c("tree_ou_index", class(idx)))
}

#' Solve the EBNM problem using an OU tree prior
#'
#' Implements an empirical Bayes normal means solver where the prior on tip
#' states is induced by a mean-zero Ornstein-Uhlenbeck process unfolding on a
#' phylogenetic tree. The root state is drawn from the stationary distribution
#' of the OU process, and two parameters are estimated by marginal maximum
#' likelihood: the reversion rate \eqn{\alpha} and stationary variance
#' \eqn{\tau^2}.
#'
#' The model is:
#' \deqn{x_\text{root} \sim N(0, \tau^2)}
#' \deqn{x_\text{child} | x_\text{parent} \sim
#'   N(e^{-\alpha \ell} x_\text{parent},\;
#'     \tau^2 (1 - e^{-2 \alpha \ell}))}
#' \deqn{y_i | x_i \sim N(x_i, s_i^2)}
#'
#' where \eqn{\ell} is the branch length connecting parent to child.
#' Marginal likelihood and posteriors are computed via an efficient
#' upward-downward message-passing algorithm implemented in C++.
#'
#' @param x A numeric vector of observations at tree tips. If named, names
#'   must match \code{tree$tip.label}; output is always in
#'   \code{tree$tip.label} order.
#'
#' @param s A scalar standard error applied to all tips, or a vector of
#'   per-tip standard errors of the same length as \code{x}. Must be
#'   strictly positive. Vectors of any other length are an error.
#'
#' @param tree A rooted phylogenetic tree of class \code{\link[ape]{phylo}}
#'   (from package \code{ape}). Must be rooted and have positive branch
#'   lengths.
#'
#' @param tree_index An optional pre-computed index returned by
#'   \code{\link{tree_ou_precomp}}. When \code{ebnm_tree_ou} is called
#'   repeatedly on the same tree (e.g., inside an iterative algorithm),
#'   building the index once and passing it here avoids rebuilding it on
#'   every call.
#'
#' @param g_init An optional \code{tree_ou} object (or an \code{ebnm} object
#'   whose \code{fitted_g} is a \code{tree_ou}) specifying an initial or
#'   fixed prior. If \code{fix_g = FALSE}, \code{g_init$alpha} and
#'   \code{g_init$stationary_var} are used as the starting point for
#'   optimization.
#'
#' @param fix_g If \code{TRUE}, fix the prior at \code{g_init} without
#'   further optimization. Requires \code{g_init} to be specified.
#'
#' @param output A character vector of values to return. Use
#'   \code{\link{ebnm_output_default}} for defaults or
#'   \code{\link{ebnm_output_all}} for all options. The
#'   \code{"posterior_sampler"} option is not supported and will be ignored.
#'
#' @param control A named list of control parameters passed to
#'   \code{\link[stats]{optim}} during optimization of \eqn{\alpha} and
#'   \eqn{\tau^2}. Supported top-level fields are \code{par_init},
#'   \code{lower}, and \code{upper}; remaining fields are passed through to
#'   \code{control = } in \code{optim}. Ignored when \code{fix_g = TRUE}.
#'
#' @return An \code{ebnm} object. See \code{\link{ebnm}} for the structure of
#'   this object. The \code{fitted_g} element is a \code{tree_ou} object with
#'   fields \code{alpha} and \code{stationary_var}. Local false sign rates
#'   (\code{lfsr}) are computed assuming a symmetric (Gaussian) posterior.
#'
#' @seealso \code{\link{ebnm}}, \code{\link{tree_ou}},
#'   \code{\link{tree_ou_precomp}}
#'
#' @export
ebnm_tree_ou <- function(
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
    stop("Package 'ape' must be installed to use ebnm_tree_ou.")
  }

  call <- match.call()
  dat <- .tbm_align_data(tree, x, s)
  y <- dat$y
  sv <- dat$s
  s2 <- sv^2
  tip.label <- dat$tip.label

  idx <- .tou_validate_tree_index(tree, tree_index)

  if (inherits(g_init, "ebnm")) {
    g_init <- g_init[["fitted_g"]]
  }
  if (!is.null(g_init) && !inherits(g_init, "tree_ou")) {
    stop("g_init must be an object of class 'tree_ou'.")
  }
  if (fix_g && is.null(g_init)) {
    stop("g_init must be provided when fix_g = TRUE.")
  }

  if (fix_g) {
    alpha_hat <- g_init$alpha
    stationary_var_hat <- g_init$stationary_var
    fit <- tree_ou_compute_cpp(
      ntip = idx$ntip,
      total_nodes = idx$total_nodes,
      root = idx$root,
      y = y,
      s2 = s2,
      edge_len_to_child = idx$edge_len_to_child,
      children = idx$children,
      internal_post = idx$internal_post,
      internal_pre = idx$internal_pre,
      alpha = alpha_hat,
      stationary_var = stationary_var_hat
    )
    loglik_val <- fit$loglik
  } else {
    bounds <- .tou_bounds(y, sv, idx, control)
    par_init <- .tou_init_par(y, sv, idx, g_init, bounds$lower, bounds$upper)

    negloglik_fn <- function(par) {
      -tree_ou_compute_cpp(
        ntip = idx$ntip,
        total_nodes = idx$total_nodes,
        root = idx$root,
        y = y,
        s2 = s2,
        edge_len_to_child = idx$edge_len_to_child,
        children = idx$children,
        internal_post = idx$internal_post,
        internal_pre = idx$internal_pre,
        alpha = par[[1]],
        stationary_var = par[[2]]
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
    loglik_val <- -opt$value

    fit <- tree_ou_compute_cpp(
      ntip = idx$ntip,
      total_nodes = idx$total_nodes,
      root = idx$root,
      y = y,
      s2 = s2,
      edge_len_to_child = idx$edge_len_to_child,
      children = idx$children,
      internal_post = idx$internal_post,
      internal_pre = idx$internal_pre,
      alpha = alpha_hat,
      stationary_var = stationary_var_hat
    )
  }

  post_mean <- stats::setNames(fit$post_tip_mean, tip.label)
  post_var  <- stats::setNames(fit$post_tip_var, tip.label)
  post_sd   <- sqrt(post_var)

  y_named <- stats::setNames(y, tip.label)
  s_named <- stats::setNames(sv, tip.label)

  retlist <- list()

  if (data_in_output(output)) {
    retlist <- add_data_to_retlist(retlist, y_named, s_named)
  }

  if (posterior_in_output(output)) {
    posterior <- list()
    if (result_in_output(output)) {
      posterior$mean  <- post_mean
      posterior$sd    <- post_sd
      posterior$mean2 <- post_mean^2 + post_var
    }
    if (lfsr_in_output(output)) {
      posterior$lfsr <- stats::pnorm(0, abs(post_mean), post_sd)
    }
    retlist <- add_posterior_to_retlist(retlist, posterior, output, y_named)
  }

  if (g_in_output(output)) {
    retlist <- add_g_to_retlist(retlist, tree_ou(alpha_hat, stationary_var_hat))
  }

  if (llik_in_output(output)) {
    df <- if (fix_g) 0L else 2L
    retlist <- add_llik_to_retlist(retlist, loglik_val, y_named, df = df)
  }

  if (sampler_in_output(output)) {
    warning("Posterior sampler is not implemented for ebnm_tree_ou.")
  }

  as_ebnm(retlist, call)
}

.tou_validate_tree_index <- function(tree, tree_index) {
  if (is.null(tree_index)) {
    return(tree_ou_precomp(tree))
  }
  if (!inherits(tree_index, "tree_ou_index") &&
      !inherits(tree_index, "tree_bm_index")) {
    stop("tree_index must be an object returned by tree_ou_precomp().")
  }
  tree_index
}

.tou_bounds <- function(y, s, idx, control) {
  min_edge <- max(min(idx$edge_len_to_child[!is.na(idx$edge_len_to_child)]), 1e-8)
  default_lower <- c(alpha = 1e-8, stationary_var = 1e-10)
  default_upper <- c(
    alpha = 100 / min_edge,
    stationary_var = max(stats::var(y), stats::median(s^2), 1.0) * 100
  )

  lower <- default_lower
  upper <- default_upper

  if (!is.null(control$lower)) {
    lower <- rep_len(control$lower, 2)
    names(lower) <- names(default_lower)
  }
  if (!is.null(control$upper)) {
    upper <- rep_len(control$upper, 2)
    names(upper) <- names(default_upper)
  }

  list(lower = lower, upper = upper)
}

.tou_init_par <- function(y, s, idx, g_init, lower, upper) {
  if (!is.null(g_init)) {
    par <- c(alpha = g_init$alpha, stationary_var = g_init$stationary_var)
  } else {
    edge_lengths <- idx$edge_len_to_child[!is.na(idx$edge_len_to_child)]
    signal_var <- max(stats::var(y) - stats::median(s^2), 1e-4)
    par <- c(
      alpha = 1 / max(stats::median(edge_lengths), 1e-6),
      stationary_var = signal_var
    )
  }

  pmin(pmax(par, lower), upper)
}
