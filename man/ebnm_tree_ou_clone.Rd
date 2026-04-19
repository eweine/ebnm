% Generated manually; update alongside R/tree_ou_clone.R
\name{ebnm_tree_ou_clone}
\alias{ebnm_tree_ou_clone}
\title{Solve the EBNM problem using a hierarchical OU forest prior}
\usage{
ebnm_tree_ou_clone(
  x,
  s = 1,
  tree,
  tree_index = NULL,
  g_init = NULL,
  fix_g = FALSE,
  output = ebnm_output_default(),
  control = NULL
)
}
\arguments{
\item{x}{A numeric vector of observations at tree tips.}

\item{s}{A scalar standard error applied to all tips, or a vector of
per-tip standard errors of the same length as \code{x}.}

\item{tree}{A rooted phylogenetic tree of class \code{\link[ape]{phylo}},
or a list of rooted \code{phylo} objects.}

\item{tree_index}{An optional pre-computed index returned by
\code{\link{tree_ou_clone_precomp}}.}

\item{g_init}{An optional \code{tree_ou_clone} object (or an \code{ebnm}
object whose \code{fitted_g} is a \code{tree_ou_clone}) specifying an
initial or fixed prior.}

\item{fix_g}{If \code{TRUE}, fix the prior at \code{g_init}.}

\item{output}{A character vector of values to return.}

\item{control}{A named list of control parameters passed to
\code{\link[stats]{optim}} during optimization.}
}
\value{
An \code{ebnm} object whose \code{fitted_g} element is a
\code{tree_ou_clone} object with fields \code{alpha},
\code{stationary_var}, and \code{mean_var}.
}
\description{
Fits an empirical Bayes normal means model on a rooted tree or forest in
which each tree has its own latent mean drawn from a mean-zero normal
distribution, and conditional on that mean each tree follows a stationary
OU process with shared parameters.
}
