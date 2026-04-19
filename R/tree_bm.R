#' Constructor for tree_bm prior class
#'
#' Creates a Brownian motion tree prior object storing the estimated diffusion
#' variance.
#'
#' @param bm_var A positive scalar: the Brownian motion variance per unit
#'   branch length.
#'
#' @return An object of class \code{tree_bm}.
#'
#' @export
#'
tree_bm <- function(bm_var) {
  if (!is.numeric(bm_var) || length(bm_var) != 1 || bm_var <= 0) {
    stop("bm_var must be a single positive number.")
  }
  structure(list(bm_var = bm_var), class = "tree_bm")
}

#' Create a flashier-compatible ebnm function for the tree BM prior
#'
#' Returns a closure with the signature \code{function(x, s, g_init, fix_g,
#' output)} expected by \pkg{flashier}. The \code{tree} and the pre-computed
#' \code{tree_index} are captured in the closure, so they are built only once
#' even when the function is called many times inside an iterative algorithm.
#'
#' @param tree A rooted phylogenetic tree of class \code{\link[ape]{phylo}},
#'   or a list of rooted \code{phylo} objects representing independent trees.
#' @param tree_index An optional pre-computed index from
#'   \code{\link{tree_bm_precomp}}. If \code{NULL} (default), it is built
#'   automatically from \code{tree}.
#'
#' @return A function suitable for use as the \code{ebnm_fn} argument to
#'   \code{\link[flashier]{flash}} or \code{\link[flashier]{flash_greedy}}.
#'
#' @seealso \code{\link{ebnm_tree_bm}}, \code{\link{tree_bm_precomp}}
#'
#' @export
#'
ebnm_tree_bm_fn <- function(tree, tree_index = NULL) {
  if (is.null(tree_index)) {
    tree_index <- tree_bm_precomp(tree)
  }
  ntip <- tree_index$ntip
  function(x, s, g_init = NULL, fix_g = FALSE, output) {
    if (length(x) != ntip) {
      # flashier validates ebnm functions by calling them with a synthetic
      # 3-element vector; fall back to a symmetric normal prior so the check
      # passes and the correct dim.signs = 0 is reported.
      return(ebnm_normal(x, s, g_init = NULL, fix_g = FALSE, output = output))
    }
    ebnm_tree_bm(
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

#' Pre-compute tree structure for repeated ebnm_tree_bm calls
#'
#' Builds the message-passing index for a phylogenetic tree or forest. Passing this
#' pre-computed index to \code{\link{ebnm_tree_bm}} via the \code{tree_index}
#' argument avoids rebuilding it on every call, which is important when
#' \code{ebnm_tree_bm} is used iteratively (e.g., as a prior in an
#' alternating optimization).
#'
#' @param tree A rooted phylogenetic tree of class \code{\link[ape]{phylo}},
#'   or a list of rooted \code{phylo} objects representing independent trees.
#'
#' @return A list of class \code{tree_bm_index} containing the pre-computed
#'   tree structure needed by \code{\link{ebnm_tree_bm}}.
#'
#' @seealso \code{\link{ebnm_tree_bm}}
#'
#' @export
#'
tree_bm_precomp <- function(tree) {
  idx <- .tree_phylo_precomp(tree)
  structure(idx, class = c("tree_bm_index", class(idx)))
}

.tree_phylo_precomp <- function(tree) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' must be installed to pre-compute tree structure.")
  }
  trees <- .tree_as_list(tree)
  tree_ntip <- vapply(trees, ape::Ntip, integer(1))
  tree_nnode <- vapply(trees, `[[`, integer(1), "Nnode")
  ntip_total <- sum(tree_ntip)
  nnode_total <- sum(tree_nnode)
  total_nodes <- ntip_total + nnode_total

  roots <- integer()
  tip_labels <- character()
  edge_len_to_child <- rep(NA_real_, total_nodes)
  children <- vector("list", total_nodes)
  internal_post <- integer()
  internal_pre <- integer()
  tip_offset <- 0L
  internal_offset <- 0L

  for (tt in seq_along(trees)) {
    tr <- trees[[tt]]
    if (!ape::is.rooted(tr)) {
      stop("Each tree must be rooted.")
    }
    if (is.null(tr$edge.length)) {
      stop("Each tree must have branch lengths.")
    }
    if (any(!is.finite(tr$edge.length)) || any(tr$edge.length <= 0)) {
      stop("All branch lengths must be finite and strictly positive.")
    }

    ntip <- tree_ntip[[tt]]
    nnode <- tree_nnode[[tt]]

    edge <- tr$edge
    edge_len <- tr$edge.length
    root <- setdiff(edge[, 1], edge[, 2])
    if (length(root) != 1L) {
      stop("Could not identify a unique root node.")
    }

    edge_len_local <- rep(NA_real_, ntip + nnode)
    children_local <- vector("list", ntip + nnode)
    for (i in seq_len(nrow(edge))) {
      p <- edge[i, 1]
      ch <- edge[i, 2]
      edge_len_local[ch] <- edge_len[i]
      children_local[[p]] <- c(children_local[[p]], ch)
    }

    tree_post <- ape::reorder.phylo(tr, order = "postorder")
    post_pars <- tree_post$edge[, 1]
    internal_post_local <- post_pars[!duplicated(post_pars, fromLast = TRUE)]

    tree_pre <- ape::reorder.phylo(tr, order = "cladewise")
    internal_pre_local <- unique(tree_pre$edge[, 1])

    map_node <- function(node_id) {
      ifelse(
        node_id <= ntip,
        tip_offset + node_id,
        ntip_total + internal_offset + (node_id - ntip)
      )
    }

    roots <- c(roots, unname(map_node(root)))
    tip_labels <- c(tip_labels, tr$tip.label)
    mapped_nodes <- unname(map_node(seq_len(ntip + nnode)))
    edge_len_to_child[mapped_nodes] <- edge_len_local
    children[mapped_nodes] <- lapply(children_local, function(ch) {
      if (length(ch) == 0L) integer() else unname(map_node(ch))
    })
    internal_post <- c(internal_post, unname(map_node(internal_post_local)))
    internal_pre <- c(internal_pre, unname(map_node(internal_pre_local)))

    tip_offset <- tip_offset + ntip
    internal_offset <- internal_offset + nnode
  }

  if (anyDuplicated(tip_labels)) {
    stop("Tip labels must be unique across all trees.")
  }

  list(
    ntip = ntip_total,
    nnode = nnode_total,
    total_nodes = total_nodes,
    root = roots,
    n_tree = length(trees),
    tip.label = tip_labels,
    edge_len_to_child = edge_len_to_child,
    children = children,
    internal_post = internal_post,
    internal_pre = internal_pre
  )
}

#' Solve the EBNM problem using a Brownian motion tree prior
#'
#' Implements an empirical Bayes normal means solver where the prior on tip
#' states is induced by a Brownian motion process unfolding on a phylogenetic
#' tree or forest. Each root state is fixed at zero, and a single variance parameter
#' \eqn{\sigma^2} (the diffusion variance per unit branch length) is estimated
#' by marginal maximum likelihood.
#'
#' The model is:
#' \deqn{x_\text{root} = 0}
#' \deqn{x_\text{child} | x_\text{parent} \sim
#'   N(x_\text{parent},\; \sigma^2 \cdot \ell)}
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
#'   (from package \code{ape}), or a list of rooted \code{phylo} objects.
#'   Trees must be rooted and have positive branch lengths. Different trees
#'   are modeled as prior-independent.
#'
#' @param tree_index An optional pre-computed index returned by
#'   \code{\link{tree_bm_precomp}}. When \code{ebnm_tree_bm} is called
#'   repeatedly on the same tree (e.g., inside an iterative algorithm),
#'   building the index once and passing it here avoids rebuilding it on
#'   every call.
#'
#' @param g_init An optional \code{tree_bm} object (or an \code{ebnm} object
#'   whose \code{fitted_g} is a \code{tree_bm}) specifying an initial or
#'   fixed prior. If \code{fix_g = FALSE}, \code{g_init$bm_var} is used as
#'   the starting point for optimization.
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
#'   \code{\link[stats]{optimize}} during optimization of \eqn{\sigma^2}.
#'   Supported fields are \code{lower}, \code{upper}, and \code{tol}.
#'   Ignored when \code{fix_g = TRUE}.
#'
#' @return An \code{ebnm} object. See \code{\link{ebnm}} for the structure of
#'   this object. The \code{fitted_g} element is a \code{tree_bm} object with
#'   a single field \code{bm_var}. Local false sign rates (\code{lfsr}) are
#'   computed assuming a symmetric (Gaussian) posterior.
#'
#' @seealso \code{\link{ebnm}}, \code{\link{tree_bm}},
#'   \code{\link{tree_bm_precomp}}
#'
#' @export
#'
ebnm_tree_bm <- function(
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
    stop("Package 'ape' must be installed to use ebnm_tree_bm.")
  }

  call <- match.call()

  # Validate and align data to tree tip order
  dat <- .tbm_align_data(tree, x, s)
  y <- dat$y
  sv <- dat$s
  s2 <- sv^2
  tip.label <- dat$tip.label

  # Use pre-computed index if supplied, otherwise build it now
  if (!is.null(tree_index)) {
    if (!inherits(tree_index, "tree_bm_index")) {
      stop("tree_index must be an object returned by tree_bm_precomp().")
    }
    idx <- tree_index
  } else {
    idx <- tree_bm_precomp(tree)
  }

  # Validate g_init
  if (inherits(g_init, "ebnm")) {
    g_init <- g_init[["fitted_g"]]
  }
  if (!is.null(g_init) && !inherits(g_init, "tree_bm")) {
    stop("g_init must be an object of class 'tree_bm'.")
  }
  if (fix_g && is.null(g_init)) {
    stop("g_init must be provided when fix_g = TRUE.")
  }

  # Determine bm_var: fixed or optimized
  if (fix_g) {
    bm_var_hat <- g_init$bm_var
    fit <- tree_bm_compute_cpp(
      ntip = idx$ntip,
      total_nodes = idx$total_nodes,
      roots = idx$root,
      y = y,
      s2 = s2,
      edge_len_to_child = idx$edge_len_to_child,
      children = idx$children,
      internal_post = idx$internal_post,
      internal_pre = idx$internal_pre,
      bm_var = bm_var_hat
    )
    loglik_val <- fit$loglik
  } else {
    loglik_fn <- function(bm_var) {
      tree_bm_compute_cpp(
        ntip = idx$ntip,
        total_nodes = idx$total_nodes,
        roots = idx$root,
        y = y,
        s2 = s2,
        edge_len_to_child = idx$edge_len_to_child,
        children = idx$children,
        internal_post = idx$internal_post,
        internal_pre = idx$internal_pre,
        bm_var = bm_var
      )$loglik
    }

    lower_bm_var <- 1e-10
    upper_bm_var <- .tbm_upper_bm_var(y, idx)

    opt_args <- c(
      list(
        f = loglik_fn,
        lower = lower_bm_var,
        upper = upper_bm_var,
        maximum = TRUE
      ),
      if (!is.null(control)) control else list()
    )
    opt <- do.call(stats::optimize, opt_args)

    bm_var_hat <- opt$maximum
    loglik_val <- opt$objective

    # Final posterior at optimum
    fit <- tree_bm_compute_cpp(
      ntip = idx$ntip,
      total_nodes = idx$total_nodes,
      roots = idx$root,
      y = y,
      s2 = s2,
      edge_len_to_child = idx$edge_len_to_child,
      children = idx$children,
      internal_post = idx$internal_post,
      internal_pre = idx$internal_pre,
      bm_var = bm_var_hat
    )
  }

  # Posterior summaries
  post_mean <- stats::setNames(fit$post_tip_mean, tip.label)
  post_var  <- stats::setNames(fit$post_tip_var, tip.label)
  post_sd   <- sqrt(post_var)

  # Named data vectors (in tip.label order) for output helpers
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
    retlist <- add_g_to_retlist(retlist, tree_bm(bm_var_hat))
  }

  if (llik_in_output(output)) {
    df <- if (fix_g) 0L else 1L
    retlist <- add_llik_to_retlist(retlist, loglik_val, y_named, df = df)
  }

  if (sampler_in_output(output)) {
    warning("Posterior sampler is not implemented for ebnm_tree_bm.")
  }

  return(as_ebnm(retlist, call))
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.tbm_align_data <- function(tree, y, s) {
  trees <- .tree_as_list(tree)
  if (any(!vapply(trees, ape::is.rooted, logical(1)))) {
    stop("Each tree must be rooted.")
  }

  tip_label <- unlist(lapply(trees, `[[`, "tip.label"), use.names = FALSE)
  if (anyDuplicated(tip_label)) {
    stop("Tip labels must be unique across all trees.")
  }

  n <- length(tip_label)
  s_scalar <- (length(s) == 1L)
  if (s_scalar) s <- rep(s, n)

  if (length(y) != n || length(s) != n) {
    stop("y and s must each have length equal to the number of tips in tree.")
  }

  if (!is.null(names(y))) {
    if (anyDuplicated(names(y))) {
      stop("names(x) must be unique.")
    }
    if (!all(tip_label %in% names(y))) {
      stop("All tree tip labels must appear in names(x).")
    }
    # Permutation that puts y into tree$tip.label order
    perm <- match(tip_label, names(y))
    y <- y[perm]
    if (!is.null(names(s)) && !s_scalar) {
      # s is a named vector: reorder it independently by name
      if (anyDuplicated(names(s))) {
        stop("names(s) must be unique.")
      }
      if (!all(tip_label %in% names(s))) {
        stop("All tree tip labels must appear in names(s).")
      }
      s <- s[tip_label]
    } else {
      # s is a scalar (already expanded) or unnamed vector: apply same
      # permutation as y so that s[i] stays paired with its y[i]
      s <- s[perm]
    }
  } else if (!is.null(names(s)) && !s_scalar) {
    stop("If s is a named vector, y must also be named.")
  } else {
    names(y) <- tip_label
    names(s) <- tip_label
  }

  if (any(!is.finite(y))) stop("x contains non-finite values.")
  if (any(!is.finite(s))) stop("s contains non-finite values.")
  if (any(s <= 0)) stop("All s must be strictly positive.")

  list(y = as.numeric(y), s = as.numeric(s), tip.label = tip_label)
}

# Upper bound for bm_var search: data variance divided by minimum edge length.
# Any bm_var beyond this would make the prior variance at every tip far exceed
# the observed variance, so the MLE cannot lie above it.
.tbm_upper_bm_var <- function(y, idx) {
  valid_len <- idx$edge_len_to_child[!is.na(idx$edge_len_to_child)]
  min_edge  <- max(min(valid_len), 1e-8)
  max(stats::var(y) / min_edge, 1.0)
}

.tree_as_list <- function(tree) {
  if (inherits(tree, "phylo")) {
    return(list(tree))
  }
  if (!is.list(tree) || length(tree) == 0L || !all(vapply(tree, inherits, logical(1), "phylo"))) {
    stop("tree must be a 'phylo' object or a non-empty list of 'phylo' objects.")
  }
  tree
}
