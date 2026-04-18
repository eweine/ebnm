# here, I just want to play around with the tip BM prior a bit
# and make sure that things are lining up with what I would expect

library(ape)

set.seed(1)
tree <- ape::rtree(75000)
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
theta_true <- setNames(x_node[seq_len(ntip)], tree$tip.label)

s_vec    <- setNames(runif(ntip, 3, 6), tree$tip.label)
y        <- setNames(rnorm(ntip, theta_true, s_vec), tree$tip.label)

fit <- ebnm_tree_bm(y, s_vec, tree)

plot(fit$posterior$mean, theta_true)

