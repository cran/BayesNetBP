
ElimTreeNodes <- function(graph, eo) {
  cluster.sets <- list()
  dag.graph <- igraph.from.graphNEL(graph)
  eo.rev <- rev(eo)
  
  # iterate over all nodes in Bayesian network
  for (i in length(eo):1){
    node <- eo.rev[i]
    neighbors <- names(neighbors(dag.graph, node, mode="all"))
    formers <- eo.rev[1:i] # nodes appear later than this node in the EO
    # a cluster for a node is formed by all its neighbors appearing later in the EO and itself
    cluster <- c(node, intersect(formers, neighbors))
    cluster.sets[[i]] <- cluster
  }
  
  # Each cluster correspond to a node in the BN, which is the eliminate node of that cluster
  names(cluster.sets) <- eo.rev
  return(cluster.sets)
}