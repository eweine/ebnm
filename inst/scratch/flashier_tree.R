library(ape)
library(ebnm)
library(flashier)

n_cells <- 1000
n_genes <- 250

set.seed(1)
tree <- ape::rtree(n_cells)
tree$edge.length <- pmax(tree$edge.length, 0.05)
ntip <- ape::Ntip(tree)

# Simulate tip states under BM with root fixed at 0
bm_var_true <- 1.0
edge_pre <- ape::reorder.phylo(tree, "cladewise")$edge
elen_pre  <- ape::reorder.phylo(tree, "cladewise")$edge.length
root_node <- setdiff(tree$edge[, 1], tree$edge[, 2])
x_node    <- rep(0, ntip + tree$Nnode)
for (i in seq_len(nrow(edge_pre))) {
  p <- edge_pre[i, 1]; ch <- edge_pre[i, 2]
  x_node[ch] <- rnorm(1, x_node[p], sqrt(bm_var_true * elen_pre[i]))
}
l_true <- setNames(x_node[seq_len(ntip)], tree$tip.label)

f_true <- rexp(n_genes)
f_true[sample(n_genes, round(n_genes / 2))] <- 0

mean_true <- tcrossprod(l_true, f_true)
E <- matrix(
  nrow = n_cells,
  ncol = n_genes,
  data = rnorm(n_cells * n_genes)
)

Y <- mean_true + E

flash_fit <- flash(
  data = Y,
  S = 1,
  ebnm_fn = list(
    ebnm::ebnm_tree_bm_fn(tree),
    ebnm::ebnm_point_exponential
  )
)

plot(flash_fit$L_pm[,1], l_true)
plot(flash_fit$F_pm[,1], f_true)



