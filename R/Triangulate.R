
Triangulate <- function(graph, eo){
  dag.graph <- igraph.from.graphNEL(graph, weight=FALSE)
  # iterate over all nodes by elimination order
  for (i in 1:length(eo)){
    node <- eo[i]
    neighbors <- names(neighbors(dag.graph, node, mode="all"))
    n.neighbors <- length(neighbors)
    if (n.neighbors>=2){
      # iterate over all pairs of neighbors
      for (p1 in 1:(n.neighbors-1)){
        for (p2 in (p1+1):n.neighbors){
          # if both of the neighbor appears later in EO & are not connected
          # then connect them
          if (which(eo==neighbors[p1])>i){
            if (which(eo==neighbors[p2])>i){
              if(!are_adjacent(dag.graph, neighbors[p1], neighbors[p2])){
                dag.graph <- add_edges(dag.graph, c(neighbors[p1],neighbors[p2]))
              }
            }
          }
        }
      }
    }
  }
  dag.tri <- igraph.to.graphNEL(dag.graph) 
  return(dag.tri)
}