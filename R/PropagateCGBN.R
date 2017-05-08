
#' Propagate the cluster tree
#' 
#' This function propagates the discrete compartment of a \code{clustertree} object.
#' 
#' @details The discrete compartment must be propagted to get the joint distributions
#' of discrete variables in each discrete clusters. A cluster tree must be propagated
#' before absorbing evidence and making queries. 
#'
#' @param tree an initialized \code{clustertree} object
#' 
#' @return a \code{clustertree} object
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
#' tree.init$propagated
#' tree.init.p <- PropagateDBN(tree.init)
#' tree.init.p$propagated
#'
#' @export

PropagateDBN <- function(tree) {
  
  post.tree <- AbsorbEvidence(tree, vars=c(), values=list())
  
}
