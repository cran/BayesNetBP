
EliminationOrder <- function(graph, discrete){
  
  dag.graph <- igraph.from.graphNEL(graph)
  
  dis.nodes <- names(discrete[discrete])
  cont.nodes <- names(discrete[!discrete])
  
  # topological order for discrete subgraph
  graph.dis <- induced_subgraph(dag.graph, dis.nodes)
  # topological order for continuous subgraph
  graph.cont <- induced_subgraph(dag.graph, cont.nodes)
  
  eo.dis <- names(topological.sort(graph.dis))
  eo.cont <- names(topological.sort(graph.cont))
  
  # eliminate continuous node first, then discrete ones
  eo <- c(rev(eo.cont), rev(eo.dis))
  return(eo)
}