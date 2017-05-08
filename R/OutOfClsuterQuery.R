

query.ooc <- function(tree, vars){
  cs.tree <- tree$cluster.tree
  tree.graph <- igraph.from.graphNEL(tree$tree.graph)
  
  cs.sets <- list()
  nm.sets <- c()
  j <- 1
  for(i in 1:length(cs.tree)) {
    if(cs.tree[[i]]@discrete) {
      cs.sets[[j]] <- cs.tree[[i]]@members
      nm.sets[j] <- cs.tree[[i]]@name
      j <- j+1
    }
  }
  
  ############
  
  vars.temp <- vars
  
  cs <- c()
  while(length(vars.temp)>0) {
    maxl <- 0
    
    inter.temp <- c()
    temp <- character(0)
    for (i in 1:length(cs.sets)) {
      
      inter <- intersect(vars.temp, cs.sets[[i]])
      if (length(inter)>=maxl) {
        maxl <- length(inter)
        temp <- nm.sets[i]
        inter.temp <- inter
      }
    }
    cs <- c(cs, temp)
    vars.temp <- setdiff(vars.temp, inter.temp)
  }
  
  # cs
  
  sub.memb <- c()
  
  for (i in 2:length(cs)) {
    path <- all_simple_paths(tree.graph, cs[1], cs[i], mode="all")[[1]]
    sub.memb <- union(sub.memb, names(path))
  }
  
  sub.graph <- induced_subgraph(tree.graph, sub.memb)
  
  # x11()
  # plot(sub.graph)
  
  node <- sub.memb[1]
  
  ooc <- list(tree=cs.tree[sub.memb], sub.graph=sub.graph, active=c(), 
              nom=cs.tree[[node]]@joint, denom=list())
  
  obj <- Distribute.OOC(ooc, node)
  temp.pot <- factor.divide(obj$nom, obj$denom)
  jnt <- marginalize.discrete(temp.pot, vars)
  
  return(jnt)
}


Distribute.OOC <- function(object.ooc, node) {
  ngbs <- neighbors(object.ooc$sub.graph, node, mode = "all")$name
  inactive <- setdiff(ngbs, object.ooc$active)
  object.ooc$active <- c(object.ooc$active, node)
  cluster <- object.ooc$tree[[node]]
  # cat(cluster@members, "\n")
  if (length(inactive)>0) {
    for (i in 1:length(inactive)) {
      # cat(node, "->", inactive[i], "\n")
      this.cluster <- object.ooc$tree[[inactive[i]]]
      object.ooc$nom <- factor.product(object.ooc$nom, this.cluster@joint, normalize=FALSE)
      # cat(names(object.ooc$nom$cpt), "\n")
      separator <- intersect(cluster@members, this.cluster@members)
      # cat(separator, "\n")
      margin <- marginalize.discrete(cluster@joint, separator)
      object.ooc$denom <- factor.product(object.ooc$denom, margin, normalize=FALSE)
      # object.ooc$joint <- factor.divide(this.cluster@joint, margin)
      object.ooc <- Distribute.OOC(object.ooc, inactive[i])
    }
  }
  return(object.ooc)
}





