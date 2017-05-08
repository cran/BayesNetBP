
SemiEliminationTree <- function(graph, cs, node.class, eo){
  
  cluster.sets <- cs
  cs.graph <- igraph.from.graphNEL(graph)
  i <- 1
  
  while (i < length(cluster.sets)){
    # cat("one iteration")
    cs.names <- names(cluster.sets)
    cluster <- cluster.sets[[i]]
    
    this.name <- cs.names[i]
    this.parent <- neighbors(cs.graph, v=this.name, mode="in")$name
    
    subset <- FALSE
    is.discrete <- prod(node.class[cluster]) # if this cluster is dicrete (all members are discete)
    if (is.discrete){ 
      # if it's discrete, search for the following clusters
      for(j in (i+1):length(cluster.sets)){
        if (prod(node.class[cluster.sets[[j]]])){ 
          # if jth cluster is discrete
          if (is.subset(cluster, cluster.sets[[j]])){ 
            # and it contains this cluster
            subset <-  TRUE 
            # set the indicator true
            break 
            # stop searching
          }
        }
      }
    }
    
    # if can not found such a cluster j, then proceed to next iteration
    if(!subset){
      i <- i+1
      next
    }
    
    # if such cluster j is found, do the following:
    
    # find all the children of this cluster in the elimination tree
    cluster.children <- neighbors(cs.graph, v=cs.names[i], mode="out")$name
    
    # find a best child among them, 
    # whose elimination node should appear as late as possible in the EO
    best.child <- 0
    best.child.name <- NULL
    minenum <- 0
    
    for(j in 1:length(cluster.children)){
      this.child <- cluster.children[[j]]
      child.mem <- cluster.sets[[this.child]]
      is.discrete <- prod(node.class[child.mem])
      if (is.discrete & is.subset(cluster, child.mem)) {
        enum <- which(eo == cs.names[i])
        
        # if the elimination node of this child appears later
        # then update the best child
        if(enum > minenum){
          
          minenum <- enum ##
          
          best.child <- j;
          best.child.name <- cluster.children[j]
        }
      }
    }
    
    # merge this cluster into its best child by
    # 1. deleting this cluster
    cs.graph <- delete_vertices(cs.graph, this.name)
    
    # 2. update parents and children of involved nodes
    if(length(this.parent)!=0){ # if this cluster has a parent
      for (j in 1:length(this.child)){
        if(this.child!=best.child.name){
          # let the parent of all its child be its best child, except for the best child itself
          cs.graph <- add_edges(cs.graph, c(best.child.name, this.child[j]))
        } else {
          # let its parent be the parent of its best child
          cs.graph <- add_edges(cs.graph, c(this.parent, this.child[j]))
        }
      }
    } else { # if it does not have a parent, then update as usual, 
      # but no need to update the parent of its best child
      for (j in 1:length(cluster.children)){
        if(cluster.children[j]!=best.child.name){
          cs.graph <- add_edges(cs.graph, c(best.child.name, cluster.children[j]))
        }
      }
    }
    
    # update the name vector of the tree
    best.child.index <- which(cs.names==best.child.name)
    cluster.sets[this.name] <- cluster.sets[best.child.index]
    cs.names[which(cs.names==this.name)] <- best.child.name
    
    cluster.sets <- cluster.sets[-best.child.index]
    cs.names <- cs.names[-best.child.index]
    names(cluster.sets) <- cs.names
    
  }
  
  # V(cs.graph)$name <- unlist(lapply(cluster.sets, paste0, collapse=","))
  # x11(); plot(igraph.to.graphNEL(cs.graph))
  
  return(igraph.to.graphNEL(cs.graph))
  
}

is.subset <- function(x,y){
  if (length(setdiff(x,y))==0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
