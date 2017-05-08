
#' Absorb evidence into the model
#'
#' @details Absorb multiple types and pieces of evidences into a \code{clustertree}
#' object. The discrete compartment of the \code{clustertree} will be automatically
#' propagated after evidence absorption, so that the object will be ready for making
#' queries and absorbing additional evidence.
#'
#' @param tree a \code{clustertree} object
#' @param vars a \code{vector} of the names of observed variables
#' @param values a \code{list} of observed values of the variables. Aside from a single value,
#' The element of the list can also be a vector of likelihood values
#' 
#' @return \code{clustertree} object with the evidence absorbed
#' 
#' @author Han Yu
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. \cr
#' \cr
#' Lauritzen, S. L., & Spiegelhalter, D. J. (1988). Local computations with probabilities on 
#' graphical structures and their application to expert systems. Journal of the Royal Statistical 
#' Society. Series B (Methodological), 157-224.
#' 
#' @import gRbase stats graphics utils
#' @importFrom igraph igraph.from.graphNEL igraph.to.graphNEL V
#' 
#' 
#' @examples 
#' 
#' data(liver)
#' cst <- ClusterTreeCompile(dag=liver$dag, node.class=liver$node.class)
#' models <- LocalModelCompile(data=liver$data, dag=liver$dag, node.class=liver$node.class)
#' tree.init <- ElimTreeInitialize(tree=cst$tree.graph, 
#'                                 dag=cst$dag, 
#'                                 model=models, 
#'                                 node.sets=cst$cluster.sets, 
#'                                 node.class=cst$node.class)
#' tree.init.p <- PropagateDBN(tree.init)
#' tree.post <- AbsorbEvidence(tree.init.p, c("Nr1i3", "chr1_42.65"), list(1,"1"))
#' 
#' @seealso \code{\link{clustertree}}
#' 
#' @export

AbsorbEvidence <- function(tree, vars, values) {
  
  ## Push continuous variables toward boundary and get discrete potentials ready
  
  ## tree <- post.tree; tree.graph <- semiET
  object <- tree
  node.class <- object$node.class
  tree <- object$cluster.tree
  
  hard <- c()
  soft <- c()
  hard.values <- list()
  soft.values <- list()
  
  if(sum(vars %in% object$absorbed.variables)!=0) {
    var.in <- vars[vars %in% object$absorbed.variables]
    msg1 <- paste0(var.in, collapse=", ")
    stop(paste0(msg1, " is/are already observed."))
  }
  
  if(sum(vars %in% object$absorbed.soft.variables)!=0) {
    var.in <- vars[vars %in% object$absorbed.variables]
    msg1 <- paste0(var.in, collapse=", ")
    warning(paste0(msg1, " has/have absorbed likelihood evidence multiple times."))
  }
  
  tree.graph <- object$tree.graph
  
  if(length(vars)!=0){
  
  var.class <- node.class[vars]
  
  for(i in 1:length(vars)) {
    if (var.class[i]) {
      if (length(values[[i]])==1){
        tree <- DiscreteEvidence(tree, vars[i], values[[i]])
        hard <- c(hard, vars[i]) #
        hard.values <- append(hard.values, values[[i]]) #
      } 
      if (length(values[[i]])>1) {
        tree <- VirtualEvidence(tree, vars[i], values[[i]])
        soft <- c(soft, vars[i]) #
        soft.values <- append(soft.values, values[i]) #
      }
      
    } 
  }
  
  for(i in 1:length(vars)) {
    if (!var.class[i]) {
      tree <- PushEvidence(tree, vars[i], values[[i]])
      hard <- c(hard, vars[i]) #
      hard.values <- append(hard.values, values[[i]]) #
    }
  }
  
  }
  
  ###################
  ## Propagate
  ###################
  
  tree.graph <- igraph.from.graphNEL(tree.graph)
  
  discrete.sets <- list()
  k <- 1
  
  discrete.clusters <- c()
  for (i in 1:length(tree)) {
    if (tree[[i]]@discrete) {
      discrete.clusters <- c(discrete.clusters, tree[[i]]@name)
      discrete.sets[[k]] <- tree[[i]]@members
      k <- k+1
    }
  }
  
  names(discrete.sets) <- discrete.clusters
  
  ## No discrete cluster
  
  if (length(discrete.clusters)==0) {
    object$cluster.tree <- tree
    object$absorbed.variables <- c(object$absorbed.variables, hard)
    object$absorbed.soft.variables <- c(object$absorbed.soft.variables, soft)
    object$absorbed.values <- append(object$absorbed.values, hard.values)
    object$absorbed.soft.values <- append(object$absorbed.soft.values, soft.values)
    names(object$absorbed.values) <- object$absorbed.variables
    names(object$absorbed.soft.values) <- object$absorbed.soft.variables
    object$propagated <- TRUE
    return(object)
  }
  
  ## single discrete cluster
  
  if (length(discrete.clusters)==1) {
    
    this.cluster <- discrete.clusters[1]
    this.tree <- tree[[this.cluster]]
    pot <- list(cpt=this.tree@cpt@table,
                prob=this.tree@cpt@logprob)
    tree[[this.cluster]]@joint <- pot
    
    object$cluster.tree <- tree
    object$absorbed.variables <- c(object$absorbed.variables, hard)
    object$absorbed.soft.variables <- c(object$absorbed.soft.variables, soft)
    object$absorbed.values <- append(object$absorbed.values, hard.values)
    object$absorbed.soft.values <- append(object$absorbed.soft.values, soft.values)
    names(object$absorbed.values) <- object$absorbed.variables
    names(object$absorbed.soft.values) <- object$absorbed.soft.variables
    object$propagated <- TRUE
    
    return(object)
  }
  
  ## Multiple discrete clusters
  
  tree.sub.graph <- induced_subgraph(tree.graph, discrete.clusters)
  potentials.sub <- list()
  for(i in 1:length(discrete.clusters)) {
    this.cluster <- discrete.clusters[i]
    this.tree <- tree[[this.cluster]]
    pot <- list(cpt=this.tree@cpt@table,
                prob=this.tree@cpt@logprob)
    potentials.sub[[i]] <- pot
  }
  names(potentials.sub) <- discrete.clusters
  
  # cat("Propagating...", "\n")
  joint.tables <- Propagate(tree.sub.graph, potentials.sub, discrete.sets)
  
  for (i in 1:length(joint.tables)) {
    this.cluster <- names(joint.tables)[i]
    tree[[this.cluster]]@joint <- joint.tables[[i]]
  }
  
  ###
  
  object$cluster.tree <- tree
  object$absorbed.variables <- c(object$absorbed.variables, hard)
  object$absorbed.soft.variables <- c(object$absorbed.soft.variables, soft)
  object$absorbed.values <- append(object$absorbed.values, hard.values)
  object$absorbed.soft.values <- append(object$absorbed.soft.values, soft.values)
  names(object$absorbed.values) <- object$absorbed.variables
  names(object$absorbed.soft.values) <- object$absorbed.soft.variables
  object$propagated <- TRUE
  
  ###
  return(object)
}

###########################################
## Reconstruct tree graph from cluster tree
###########################################

clustertree.to.graph <- function(tree) {
  
  from <- c()
  to <- c()
  
  for (i in 1:length(tree)) {
    
    from <- c(from, rep(tree[[i]]@name, length(tree[[i]]@children)))
    to <- c(to, tree[[i]]@children)
    edge.list <- data.frame(from, to)
    graph <- graph_from_data_frame(edge.list, directed=TRUE)
    
  }
  
  tree.dag <- igraph.to.graphNEL(graph)
  return(tree.dag)
  
}


###########################################
## Propagate
###########################################

Propagate <- function(tree.graph, potentials, cluster.sets){
  
  # tree.graph <- tree.sub.graph; potentials <- potentials.sub; cluster.sets <- discrete.sets
  
  cluster.tree <- list(
    # bn=dag,
    tree=tree.graph,
    clusters=cluster.sets, 
    # assignment=asgn,
    collected=c(), active=c(), potentials=potentials, joint=potentials)
  
  clusters <- names(potentials)
  result <- list()
  
  ## NEW version of getting joints
  
  # collect
  ce <- CollectEvidence(cluster.tree, clusters[1])
  # reset active nodes
  ce$active <- c()
  # distribute
  de <- DistributeEvidence(ce, clusters[1])
  result <- de$joint
  
  ## old version of getting joints ###################
  if(FALSE){
    for (i in 1:length(clusters)){
      ce <- CollectEvidence(cluster.tree, clusters[i])
      result[[i]] <- ce$potentials[[clusters[i]]]
      names(result) <- clusters
    }
  }
  ####################################################
  return(result)
}

# ppgt <- Propagate(tree.sub.graph, potentials.sub, discrete.sets)

###################################################

Absorb <- function(absorbedTo, absorbedFrom, separator, distribute=FALSE){
  pot1 <- absorbedFrom
  pot2 <- absorbedTo
  
  # pot2 <- cluster.tree$potentials[["Cyp2b10"]]; pot1 <- cluster.tree$potentials[["HDL"]];
  # separator <- intersect( cluster.tree$clusters[["Cyp2b10"]],  cluster.tree$clusters[["HDL"]])
  
  inter.var <- separator
  sep <- marginalize.discrete(pot1, inter.var)
  
  results <- list()
  if (distribute) {
    results[[1]] <- NULL
  } else {
    results[[1]] <- factor.divide(pot1, sep)
  }
  
  results[[2]] <- factor.product(pot2, sep)
  return(results)
}

#Absorb(cluster.tree$potentials[["Cyp2b10"]], cluster.tree$potentials[["HDL"]],
#      separator = intersect( cluster.tree$clusters[["Cyp2b10"]],  cluster.tree$clusters[["HDL"]]))


###########################################
## Collect evidence
###########################################

CollectEvidence <- function(cluster.tree, node) {
  clique.names <- names(V(cluster.tree$tree))
  ngbs <- neighbors(cluster.tree$tree, node, mode = "all")$name
  inactive <- setdiff(ngbs, cluster.tree$active)
  cluster.tree$active <- c(cluster.tree$active, node)
  
  for (ngb in inactive){
    # cat("collecting for ", node, "from", ngb)
    cluster.tree <- CollectEvidence(cluster.tree, ngb)
    # collected <- ce[[1]]
    # cst <- ce[[2]]
  }
  
  if (length(inactive)>0) {
    for (i in 1:length(inactive)) {
      abb <- Absorb(cluster.tree$potentials[[node]], cluster.tree$potentials[[inactive[i]]],
                    separator = intersect( cluster.tree$clusters[[node]],  cluster.tree$clusters[[inactive[i]]]))
      cluster.tree$potentials[[node]] <- abb[[2]]
      cluster.tree$potentials[[inactive[i]]] <- abb[[1]]
      # cat(inactive[i], " -> ", node, "\n")
    }
    
  }
  cluster.tree$collected <- c(cluster.tree$collected, node)
  # customizedPlot(graph, collect.update)
  return(cluster.tree)
}

###########################################
## Distribute evidence
###########################################

# need to reset the active nodes of cluster.tree after collecting evidence

DistributeEvidence <- function(cluster.tree, node){
  clique.names <- names(V(cluster.tree$tree))
  ngbs <- neighbors(cluster.tree$tree, node, mode = "all")$name
  
  inactive <- setdiff(ngbs, cluster.tree$active)
  cluster.tree$active <- c(cluster.tree$active, node)
  
  cluster.tree$joint[[node]] <- cluster.tree$potentials[[node]]
  
  if (length(inactive)>0) {
    for (i in 1:length(inactive)) {
      abb <- Absorb(cluster.tree$potentials[[inactive[i]]], cluster.tree$potentials[[node]],
                    separator = intersect( cluster.tree$clusters[[node]],  cluster.tree$clusters[[inactive[i]]]),
                    distribute = TRUE)
      # cluster.tree$potentials[[node]] <- abb[[1]]
      cluster.tree$potentials[[inactive[i]]] <- abb[[2]]
      cluster.tree <- DistributeEvidence(cluster.tree, inactive[i])
    }
  }
  
  return(cluster.tree)
}

###########################################
## Add Evidence to LPPotential
###########################################

addEvidence <- function(theObject, vars, vals)
{
  varinds <- c()
  for(i in 1:length(vars)) {
    varinds[i] <- which(theObject@tail==vars[i])
  }
  
  if (length(varinds)==0) {
    return(theObject)
  }
  
  theObject@const <- theObject@const + 
    sum(theObject@params[varinds] * vals)
  theObject@tail <- theObject@tail[-varinds]
  theObject@params <- theObject@params[-varinds]
  return(theObject)
}

###########################################
## Absorb Continuous Evidence
###########################################

PushEvidence <- function(tree.push, var, val){
  # Step 1
  # tree.push <- tree.init; var <- "Neu1"; val <- 10
  for (i in length(tree.push):1){
    if (tree.push[[i]]@name == var){
      break
    }
    if (tree.push[[i]]@activeflag){
      for (j in 1:length(tree.push[[i]]@lppotential)){
        
        if(var %in% tree.push[[i]]@lppotential[[j]]@tail){
          tree.push[[i]]@lppotential[[j]] <- 
            addEvidence(tree.push[[i]]@lppotential[[j]], var, val)
          
        }
        # print(tree.push[[i]]@lppotential[[j]]@tail)
      }
    }
  }
  
  # Step 2
  n <- i
  tree.push[[n]]@postbag <- tree.push[[n]]@lppotential
  tree.push[[n]]@lppotential <- list()
  tree.push[[n]]@activeflag <- FALSE
  
  # cat(paste0(rep("*",50), collapse=""), "\n")
  # cat("Tree with initialized pushing loop", "\n")
  # cat(paste0(rep("*",50), collapse=""), "\n")
  # PrintTree(tree.push)
  
  # Step 3
  
  
  while((length(tree.push[[n]]@parent)!=0) & 
        !tree.push[[tree.push[[n]]@parent]]@discrete){
    tree.push[[tree.push[[n]]@parent]]@postbag <- tree.push[[n]]@postbag
    
    #########################################
    flag <- tree.push[[tree.push[[n]]@parent]]@activeflag
    if (length(tree.push[[tree.push[[n]]@parent]]@lppotential)==0){
      flag <- FALSE
    } else {
      lp.head <- tree.push[[tree.push[[n]]@parent]]@lppotential[[1]]@head
      postbag.tail <- tree.push[[tree.push[[n]]@parent]]@postbag[[1]]@tail
      if( !lp.head %in% postbag.tail){
        flag <- FALSE
      }
    }
    #########################################
    
    #if (tree.push[[tree.push[[n]]@parent]]@activeflag){
    if (flag) {
      
      newBag <- BagExchange.2(tree.push[[tree.push[[n]]@parent]]@postbag,
                            tree.push[[tree.push[[n]]@parent]]@lppotential)
      tree.push[[tree.push[[n]]@parent]]@postbag <- newBag[[2]]
      tree.push[[tree.push[[n]]@parent]]@lppotential <- newBag[[1]]
      
      for (i in 1:length(tree.push[[tree.push[[n]]@parent]]@lppotential)){
        tree.push[[tree.push[[n]]@parent]]@lppotential[[i]] <- 
          addEvidence(tree.push[[tree.push[[n]]@parent]]@lppotential[[i]], var, val)
      }
    }
    
    tree.push[[n]]@postbag <- list()
    
    # cat(paste0(rep("*",50), collapse=""), "\n")
    # cat("Push Evidence: ", tree.push[[n]]@name)
    n <- which(names(tree.push) == tree.push[[n]]@parent)
    
    # cat(" ->", tree.push[[n]]@name, "\n")
    # cat(paste0(rep("*",50), collapse=""), "\n")
    # PrintTree(tree.push)
  }
  
  # step 4 v2
  
  weight <- c()
  
  vars <- tree.push[[n]]@postbag[[1]]@conditionvals
  values.vec <- c()
  for(i in 1:length(tree.push[[n]]@postbag)) {
    this.pot <- tree.push[[n]]@postbag[[i]]
    this.weight <- dnorm(val, mean = this.pot@const, sd = sqrt(this.pot@sigma))
    this.value <- this.pot@conditionvalues
    values.vec <- c(values.vec, this.value)
    weight <- c(weight, this.weight)
    
  }
  
  this.cpt <- data.frame(matrix(values.vec, ncol=length(vars), byrow=TRUE))
  colnames(this.cpt) <- vars
  
  pot <- list(cpt=this.cpt, prob=weight)
  pot.parent <- list(cpt=tree.push[[tree.push[[n]]@parent]]@cpt@table,
                     prob=tree.push[[tree.push[[n]]@parent]]@cpt@logprob)
  pot.post <- factor.product(pot, pot.parent)
  
  tree.push[[tree.push[[n]]@parent]]@cpt@factors <- names(pot.post$cpt)
  tree.push[[tree.push[[n]]@parent]]@cpt@table <- pot.post$cpt
  tree.push[[tree.push[[n]]@parent]]@cpt@logprob <- pot.post$prob
  
  
  #####################
  if(FALSE){
    
  # Step 4, old version
  weight <- c()
  for(i in 1:length(tree.push[[n]]@postbag)) {
    this.pot <- tree.push[[n]]@postbag[[i]]
    this.weight <- dnorm(val, mean = this.pot@const, sd = sqrt(this.pot@sigma))
    this.vars <- this.pot@conditionvals
    parent.cpt <- tree.push[[tree.push[[n]]@parent]]@cpt@table
    parent.prob <- tree.push[[tree.push[[n]]@parent]]@cpt@logprob
    for (j in 1:nrow(parent.cpt)) {
      this.config <- parent.cpt[j,]
      # is.compatible <- prod(this.config[this.vars]==this.pot@conditionvalues)
      
      if(length(this.vars)==0) {
        is.compatible <- TRUE
      }else{
        is.compatible <- prod(this.config[this.vars]==this.pot@conditionvalues)
      }
      
      if(is.compatible){
        weight[j] <- this.weight
      }
    }
  }
  
  tree.push[[tree.push[[n]]@parent]]@cpt@logprob <- 
    weight*tree.push[[tree.push[[n]]@parent]]@cpt@logprob
  
  tree.push[[tree.push[[n]]@parent]]@cpt@logprob <- 
    tree.push[[tree.push[[n]]@parent]]@cpt@logprob/(sum(tree.push[[tree.push[[n]]@parent]]@cpt@logprob))
  
  }
  #####################
  
  return(tree.push)
}

###########################################
## Absorb Discrete Evidence
###########################################

# tree.push <- tree;  var <- "asia"; val <- "yes";

DiscreteEvidence <- function(tree.push, var, val) {
  for (i in 1:length(tree.push)) {
    if(tree.push[[i]]@discrete) {
      this.cpt <- tree.push[[i]]@cpt
      
      if(var %in% colnames(this.cpt@table)){
        keep <- which(this.cpt@table[[var]] == val)
        this.cpt@table <- this.cpt@table[,-which(colnames(this.cpt@table)==var), drop=FALSE]
        this.cpt@table <- this.cpt@table[keep, , drop=FALSE]
        this.cpt@logprob <- this.cpt@logprob[keep]
        this.cpt@logprob <- this.cpt@logprob/sum(this.cpt@logprob)
        tree.push[[i]]@cpt <- this.cpt
      }
      
    } else {
      ## update LPPotentials
      this.lps <- tree.push[[i]]@lppotential
      temp.lps <- list()
      
      if (length(this.lps)>0) {
        for (j in 1:length(this.lps)) {
          this.lp <- this.lps[[j]]
          names(this.lp@conditionvalues) <- this.lp@conditionvals
          
          if (!var %in% this.lp@conditionvals){
            temp.lps <- c(temp.lps, this.lp)
            next
          }
          
          if (this.lp@conditionvalues[[var]] == val){
            rm.ind <- which(this.lp@conditionvals==var)
            this.lp@conditionvals <- this.lp@conditionvals[-rm.ind]
            this.lp@conditionvalues <- this.lp@conditionvalues[-rm.ind]
            temp.lps <- c(temp.lps, this.lp)
          }
        }
        tree.push[[i]]@lppotential <- temp.lps
      }
      ## update Postbag
      
      this.postbag <- tree.push[[i]]@postbag
      temp.postbag <- list()
      
      if (length(this.postbag)>0) {
        for (j in 1:length(this.postbag)) {
          this.lp <- this.postbag[[j]]
          names(this.lp@conditionvalues) <- this.lp@conditionvals
          
          if (!var %in% this.lp@conditionvals){
            temp.postbag <- c(temp.postbag, this.lp)
            next
          }
          
          if (this.lp@conditionvalues[[var]] == val){
            rm.ind <- which(this.lp@conditionvals==var)
            this.lp@conditionvals <- this.lp@conditionvals[-rm.ind]
            this.lp@conditionvalues <- this.lp@conditionvalues[-rm.ind]
            temp.postbag <- c(temp.postbag, this.lp)
          }
        }
        tree.push[[i]]@postbag <- temp.postbag
      }
    }
  }
  return(tree.push)
}


###########################################
## Virtual Evidence
###########################################

# tree <- tree.init; var <- vars[1]; val <- values[[1]];

VirtualEvidence <- function(tree, var, val) {
  for (i in 1:length(tree)){
    this.tree <- tree[[i]]
    this.table <- this.tree@cpt@table
    if (var %in% this.tree@members) {
      pot1 <- list(cpt = this.table,
                   prob = this.tree@cpt@logprob)
      
      df.temp <- data.frame(names(val))
      names(df.temp) <- var
      
      pot2 <- list(cpt = df.temp,
                   prob = val/sum(val))
      
      pot.new <- factor.product(pot1, pot2)
      
      this.tree@cpt@table <- pot.new$cpt
      this.tree@cpt@logprob <- pot.new$prob
      tree[[i]] <- this.tree
      
      break
    }
  }
  return(tree)
}


###########################################
## Exchange operation
###########################################

Exchange <- function(pot1, pot2) {
  exc.z <- pot1@head
  exc.y <- pot2@head
  
  pot1.tail <- pot1@tail
  pot2.tail <- pot2@tail
  
  names(pot1@params) <- pot1@tail
  names(pot2@params) <- pot2@tail
  
  if (!exc.y %in% pot1@tail) {
    pot1@tail <- c(pot1@tail, exc.y)
    pot1@params <- c(pot1@params, 0)
    names(pot1@params) <- pot1@tail
  }
  
  pot1.tail <- pot1@tail
  pot2.tail <- pot2@tail
  
  b <- pot1@params[exc.y]
  w <- union(pot1.tail, pot2.tail)
  w <- w[-which(w == exc.y)]
  
  if(length(w)==0){
    a <- c <- 0
  } else {
    a <- c <- rep(0, length(w))
    names(a) <- names(c) <- w
    
    
    for (i in 1:length(w)) {
      if (w[i] %in% pot1.tail) {
        a[i] <- pot1@params[which(pot1.tail==w[i])]
      }
      if (w[i] %in% pot2.tail) {
        c[i] <- pot2@params[which(pot2.tail==w[i])]
      }
    }
  }
  
  a0 <- pot1@const
  c0 <- pot2@const
  
  newpot1 <- pot1
  newpot2 <- pot2
  
  newpot2@params <- a + b*c 
  newpot2@const <- a0 + b*c0
  newpot2@sigma <- pot1@sigma + b^2*pot2@sigma
  
  den <- newpot2@sigma
  newpot1@params <- (c*pot1@sigma - a*b*pot2@sigma)/den
  newpot1@const <- (c0*pot1@sigma - a0*b*pot2@sigma)/den
  
  param.z <- b*pot2@sigma/den
  newpot1@params <- c(param.z, newpot1@params)
  newpot1@tail <- c(exc.z, w)
  names(newpot1@params) <- newpot1@tail
  newpot1@sigma <- pot1@sigma*pot2@sigma/den
  
  names(newpot1@sigma) <- NULL
  names(newpot2@sigma) <- NULL
  names(newpot1@const) <- NULL
  names(newpot2@const) <- NULL
  
  newpot1@head <- exc.y
  newpot2@head <- exc.z
  
  keep.1 <- which(newpot1@params!=0)
  keep.2 <- which(newpot2@params!=0)
  
  newpot1@params <- newpot1@params[keep.1]
  newpot2@params <- newpot2@params[keep.2]
  
  newpot1@tail <- names(newpot1@params)
  newpot2@tail <- names(newpot2@params)
  
  cond.1 <- pot1@conditionvals
  cond.2 <- pot2@conditionvals
  condv.1 <- pot1@conditionvalues
  condv.2 <- pot2@conditionvalues
  names(condv.1) <- cond.1
  names(condv.2) <- cond.2
  
  ## processing conditional variables and their values
  
  # condition 1
  
  if (length(cond.1)==0) {
    
    newpot1@conditionvals <- cond.2
    newpot2@conditionvals <- cond.2
    newpot1@conditionvalues <- condv.2
    newpot2@conditionvalues <- condv.2
    
    result <- list(newpot1, newpot2)
    return(result)
  }
  
  # condition 2
  
  if (length(cond.2)==0) {
    
    newpot1@conditionvals <- cond.1
    newpot2@conditionvals <- cond.1
    newpot1@conditionvalues <- condv.1
    newpot2@conditionvalues <- condv.1
    
    result <- list(newpot1, newpot2)
    return(result)
  }
  
  # condition 3
  
  union.cond <- union(cond.1, cond.2)
  common.val <- rep(NA, length(union.cond))
  names(common.val) <- union.cond
  
  for (i in 1:length(union.cond)){
    if(union.cond[i] %in% cond.1){
      common.val[union.cond[i]] <- condv.1[union.cond[i]]
    } else {
      common.val[union.cond[i]] <- condv.2[union.cond[i]]
    }
  }
  
  newpot1@conditionvals <- union.cond
  newpot2@conditionvals <- union.cond
  newpot1@conditionvalues <- common.val
  newpot2@conditionvalues <- common.val
  
  result <- list(newpot1, newpot2)
  return(result)
}

###########################################
## Check if two potentials are exchangeable
###########################################

is.Exchangeable <- function(pot1, pot2){
  var.1 <- pot1@conditionvals
  var.2 <- pot2@conditionvals
  var.inter <- intersect(var.1, var.2)
  if (length(var.inter)==0){
    return(TRUE)
  }
  names(pot1@conditionvalues) <- var.1
  names(pot2@conditionvalues) <- var.2
  val.1 <- pot1@conditionvalues[var.inter]
  val.2 <- pot2@conditionvalues[var.inter]
  if (prod(val.1==val.2)){
    return(TRUE) 
  } else {
    return(FALSE)
  }
}

###########################################
## Exhange two bags of potentials
###########################################

BagExchange <- function(bag1, bag2){
  output <- newbag1 <- newbag2 <- list()
  for (i in 1:length(bag1)) {
    for (j in 1:length(bag2)) {
      lp1 <- bag1[[i]]
      lp2 <- bag2[[j]]
      if (is.Exchangeable(lp1, lp2)) {
        newlps <- Exchange(lp1, lp2)
        newbag1 <- c(newbag1, newlps[[1]])
        newbag2 <- c(newbag2, newlps[[2]])
      }
    }
  }
  output[[1]] <- newbag1
  output[[2]] <- newbag2
  return(output)
}
