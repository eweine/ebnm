% Generated manually; update alongside R/tree_ou_clone.R
\name{tree_ou_clone}
\alias{tree_ou_clone}
\title{Constructor for tree_ou_clone prior class}
\usage{
tree_ou_clone(alpha, stationary_var, mean_var)
}
\arguments{
\item{alpha}{A positive scalar giving the OU reversion rate per unit
branch length.}

\item{stationary_var}{A positive scalar giving the stationary variance of
the OU process around each tree-specific mean.}

\item{mean_var}{A positive scalar giving the variance of the
tree-specific means.}
}
\value{
An object of class \code{tree_ou_clone}.
}
\description{
Creates a hierarchical OU forest prior with shared OU parameters and
tree-specific means drawn from a mean-zero normal distribution.
}
