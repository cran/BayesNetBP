#####################################
## function for initializing continuous variable 
#####################################

InitCont <- function(clustertree, node, dag, lp.bags){
  this.node <- node
  cs.tree <- clustertree
  
  
  dag.graph <- igraph.from.graphNEL(dag)
  parent.nodes <- neighbors(dag.graph, v=this.node, mode="in")$name # parents(this.node, dag)
  
  bagcell <- lp.bags[[this.node]]
  
  for(j in 1:length(cs.tree)){
    this.cst <- cs.tree[[j]]
    if (is.subset(c(this.node, parent.nodes), this.cst@members)) {
      elim.node <- this.cst@name
      break
    }
  }
  
  if(elim.node == this.node){
    cs.tree[[j]]@lppotential <- bagcell
  } else {
    cs.tree[[j]]@postbag <- bagcell
  }
  return(cs.tree)
}


#####################################
## function for assigning universe 
#####################################

assignUniverse <- function(dag, universes, nodes){
  # universes <- cs.2
  assignment <- list()
  node.names <- dag@nodes
  assigned <- c()
  assigned <- setdiff(dag@nodes, nodes)
  dag.graph <- igraph.from.graphNEL(dag)
  
  i <- 1
  for (universe in universes){
    temp <- c()
    for (node in universe){
      if (node %in% assigned) next
      node.parents <- names(neighbors(dag.graph, node, mode="in"))
      if (length(node.parents)==0) {
        temp <- c(temp, node)
        assigned <- c(assigned, node)
        next
      }
      if (prod(node.parents %in% universe) == 1) {
        temp <- c(temp, node)
        assigned <- c(assigned, node)
      }
    }
    assignment[[i]] <- temp
    i <- i+1
  }
  names(assignment) <- names(universes)
  return(assignment)
}


####################
## Factor product
####################

factor.product <- function(pot1, pot2, normalize=TRUE){
  
  p <- c()
  r <- 1
  
  cpt10 <- pot1$cpt
  cpt20 <- pot2$cpt
  
  vars1 <- names(cpt10)
  vars2 <- names(cpt20)
  inter.vars <- intersect(vars1, vars2)
  
  prob1 <- pot1$prob
  prob2 <- pot2$prob
  
  if (length(vars1)==0) { return(pot2) }
  if (length(vars2)==0) { return(pot1) }
  
  cpt1 <- cpt10[inter.vars]
  cpt2 <- cpt20[inter.vars]
  
  all.vars <- c(vars1, vars2)
  index <- match(unique(all.vars), all.vars)
  
  cpt1.char <- apply(cpt1, 1, paste0, collapse="%")
  cpt2.char <- apply(cpt2, 1, paste0, collapse="%")
  uni.char <- unique(cpt1.char)
  
  cpt.v <- c()
  
  
  for (i in 1:length(uni.char)) {
    this.char <- uni.char[i]
    pos1 <- which(cpt1.char==this.char)
    pos2 <- which(cpt2.char==this.char)
    
    for (j in 1:length(pos1)){
      for (k in 1:length(pos2)){
        
        p1 <- pos1[j]
        p2 <- pos2[k]
        
        all.vals <- c( as.character(t( cpt10[p1,] )), as.character(t( cpt20[p2,]) ))
        vals <- as.character(all.vals[index])
        cpt.v <- c(cpt.v, vals)
        
        p[r] <- prob1[p1]*prob2[p2]
        r <- r+1
        
      }
    }
  }
  cpt <- matrix(cpt.v, ncol=length(index), byrow=TRUE)
  cpt <- data.frame(cpt)
  colnames(cpt) <- all.vars[index]
  
  if (normalize) {
    p <- p/sum(p)
  }
  
  
  result <- list(cpt=cpt, prob=p)
  return(result)
}


###########################
## factor divide
###########################

factor.divide <- function(pot1, pot2){
  
  p <- c()
  r <- 1
  
  cpt10 <- pot1$cpt
  cpt20 <- pot2$cpt
  
  vars1 <- names(cpt10)
  vars2 <- names(cpt20)
  inter.vars <- intersect(vars1, vars2)
  
  prob1 <- pot1$prob
  prob2 <- pot2$prob
  
  cpt1 <- cpt10[inter.vars]
  cpt2 <- cpt20[inter.vars]
  
  all.vars <- c(vars1, vars2)
  index <- match(unique(all.vars), all.vars)
  
  cpt1.char <- apply(cpt1, 1, paste0, collapse="%")
  cpt2.char <- apply(cpt2, 1, paste0, collapse="%")
  uni.char <- unique(cpt1.char)
  
  cpt.v <- c()
  
  
  for (i in 1:length(uni.char)) {
    this.char <- uni.char[i]
    pos1 <- which(cpt1.char==this.char)
    pos2 <- which(cpt2.char==this.char)
    
    for (j in 1:length(pos1)){
      for (k in 1:length(pos2)){
        
        p1 <- pos1[j]
        p2 <- pos2[k]
        
        # all.vals <- as.character(c(t(cpt10[p1,]), t(cpt20[p2,])))
        all.vals <- c( as.character(t( cpt10[p1,] )), as.character(t( cpt20[p2,]) ))
        vals <- as.character(all.vals[index])
        cpt.v <- c(cpt.v, vals)
        
        if(prob2[p2]==0) {
          p[r] <- 0
        } else {
          p[r] <- prob1[p1]/prob2[p2]
        }
        
        r <- r+1
        
      }
    }
  }
  cpt <- matrix(cpt.v, ncol=length(index), byrow=TRUE)
  cpt <- data.frame(cpt)
  colnames(cpt) <- all.vars[index]
  
  # p <- p/sum(p)
  
  result <- list(cpt=cpt, prob=p)
  return(result)
}


########################################
## Conditional dsitribution 
## vars: the variables conditioned on
########################################

conditional <- function(pot, vars) {
  pot2 <- marginalize.discrete(pot, vars)
  pot3 <- factor.divide(pot, pot2)
  return (pot3)
}

####################
## Marginalize
####################

marginalize.discrete <- function(pot, vars) {
  
  # pot <- pot1; vars <- c("asia", "tub");
  cpt0 <- pot$cpt
  prob0 <- pot$prob
  
  pot.vars <- names(cpt0) ######
  vars <- intersect(pot.vars, vars) ######
  
  cpt.sub <- cpt0[vars]
  cpt.char <- apply(cpt.sub, 1, paste0, collapse="%")
  
  index <- match(unique(cpt.char), cpt.char)
  cpt <- data.frame(cpt.sub[index,])
  colnames(cpt) <- vars
  
  prob <- c()
  
  for (i in 1:length(unique(cpt.char))) {
    pos <- which(cpt.char==unique(cpt.char)[i])
    prob[i] <- sum(prob0[pos])
  }
  
  result <- list(cpt=cpt, prob=prob)
  return(result)
}

####################
## a helper function for checking if x is a subset of y #####
####################
is.subset <- function(x,y){
  if (length(setdiff(x,y))==0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

########################################
## Extract info from qtl fit results
########################################
#' @importFrom graph nodes
#' @importFrom qtlnet loci.qtlnet
#' @importFrom qtl scanone
extractQTL <- function(qtl.fit) {
  
  qtl.graph <- igraph.qtlnet(qtl.fit)
  qtl <- qtl.fit$cross
  
  dag <- igraph.to.graphNEL(qtl.graph)
  
  if (!is.DAG(dag)) stop("Graph is not a DAG.")
  
  graph::nodes(dag) <- gsub("@", "_", graph::nodes(dag))
  node.names <- graph::nodes(dag)
  
  # pheno <- qtl$pheno[,pheno.cols]
  pheno <- qtl$pheno
  
  loci <- qtlnet::loci.qtlnet(qtl.fit)
  locus <- unique(unlist(loci))
  
  locus <- gsub("@", "_", locus)
  qtl.df <- loci.loc(locus)
  
  discrete.nodes <- locus
  continuous.nodes <- names(pheno)
  
  node.class <- node.names %in% discrete.nodes
  names(node.class) <- node.names
  
  geno.list <- list()
  markers <- c()
  
  for(i in 1:nrow(qtl.df)) {
    markers[i] <- find.marker(qtl, qtl.df$chr[i], qtl.df$location[i])
    geno.list[[i]] <- 
      data.frame(qtl$geno[[qtl.df$chr[i]]]$data)[[markers[i]]]
  }
  
  geno <- matrix(unlist(geno.list), byrow=FALSE, ncol=length(geno.list))
  
  colnames(geno) <- locus
  
  dat <- data.frame(geno, pheno)
  
  result <- list(data=dat,
                 dag=dag,
                 node.class=node.class)
  return(result)
}

