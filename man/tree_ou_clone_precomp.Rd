% Generated manually; update alongside R/tree_ou_clone.R
\name{tree_ou_clone_precomp}
\alias{tree_ou_clone_precomp}
\title{Pre-compute tree structure for repeated ebnm_tree_ou_clone calls}
\usage{
tree_ou_clone_precomp(tree)
}
\arguments{
\item{tree}{A rooted phylogenetic tree of class \code{\link[ape]{phylo}},
or a list of rooted \code{phylo} objects representing independent trees.}
}
\value{
A list of class \code{tree_ou_clone_index}.
}
\description{
Builds the message-passing index for a phylogenetic tree or forest.
}
