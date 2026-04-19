% Generated manually; update alongside R/tree_ou_clone.R
\name{ebnm_tree_ou_clone_fn}
\alias{ebnm_tree_ou_clone_fn}
\title{Create a flashier-compatible ebnm function for the tree OU clone prior}
\usage{
ebnm_tree_ou_clone_fn(tree, tree_index = NULL)
}
\arguments{
\item{tree}{A rooted phylogenetic tree of class \code{\link[ape]{phylo}},
or a list of rooted \code{phylo} objects representing independent trees.}

\item{tree_index}{An optional pre-computed index from
\code{\link{tree_ou_clone_precomp}}.}
}
\value{
A function suitable for use as the \code{ebnm_fn} argument to
\code{\link[flashier]{flash}} or \code{\link[flashier]{flash_greedy}}.
}
\description{
Returns a closure with the signature
\code{function(x, s, g_init, fix_g, output)} expected by \pkg{flashier}.
}
