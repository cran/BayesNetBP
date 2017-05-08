#' Initialize the elimination tree
#'
#' Initialize the elimination tree with the local models
#'
#' @details Initialize the elimination tree with the local models
#'
#' @param tree a \code{graphNEL} object of the elimination tree
#' @param dag a \code{graphNEL} object of the Bayesian network
#' @param model a \code{list} of local models built from \code{LocalModelCompile} function
#' @param node.sets a \code{list} of cluster sets obtained from \code{ClusterTreeCompile} function
#' @param node.class a named \code{vector} of \code{logical} values, \code{TRUE} if node 
#' is discrete, \code{FASLE} if otherwise
#' 
#' @return \code{clustertree} object with the local models incorporated
#' 
#' @author Han Yu
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. 
#' 
#' @import qtlnet doBy
#' @importFrom graph nodes
#' @importFrom igraph neighbors
#' @importFrom methods new
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
#' 
#' @seealso The functions \code{\link{ClusterTreeCompile}} and \code{\link{LocalModelCompile}} provide necessary
#' objects to obtain \code{\link{clustertree}} object by initializing the elimination tree through this function.
#' 
#' @export

ElimTreeInitialize <- function(tree, dag, model, node.sets, node.class){
  
  # tree <- semiET; model <- models; node.sets <- cs; node.class <- node.class
  
  tree.nodes <- graph::nodes(tree)
  tree.graph <- igraph.from.graphNEL(tree)
  
  elim.tree <- list() 
  discrete.clusters <- c()
  
  discrete.nodes <- names(node.class)[node.class]
  
  # initialize with topological information
  
  for (i in 1:length(tree.nodes)) {
    this.tree <- new("ClusterSetTree", 
                     name=tree.nodes[i],
                     members = node.sets[[tree.nodes[i]]], 
                     index=i)
    this.tree@parent <- neighbors(tree.graph, v=this.tree@name, mode="in")$name
    this.tree@children <- neighbors(tree.graph, v=this.tree@name, mode="out")$name
    this.tree@activeflag <- TRUE
    this.tree@discrete <- as.logical(prod(node.class[this.tree@members]))
    
    if(this.tree@discrete) {
      discrete.clusters <- c(discrete.clusters, this.tree@name)
    }
    
    elim.tree[[i]] <- this.tree
  }
  names(elim.tree) <- tree.nodes
  
  # get discrete subgraph of the cluster tree
  
  tree.discrete <- induced_subgraph(tree.graph, discrete.clusters)
  asgn <- assignUniverse(dag, node.sets[discrete.clusters], discrete.nodes)
  
  # initialize with local models
  
  for (i in 1:length(discrete.clusters)) {
    this.cluster <- discrete.clusters[i]
    for (j in 1:length(asgn[[this.cluster]])) {
      this.asgn <- asgn[[this.cluster]][j]
      if (j==1) {
        pot <- model$pots[[this.asgn]]
      } else {
        pot <- factor.product(pot, model$pots[[this.asgn]])
      }
    }
    
    cpt <- new("CondProbTable",
               factors = elim.tree[[this.cluster]]@members,
               table = pot$cpt,
               logprob = pot$prob)
    elim.tree[[this.cluster]]@cpt <- cpt
  }
  
  for (i in 1:length(tree.nodes)) {
    this.tree <- tree.nodes[i]
    if(elim.tree[[this.tree]]@discrete){
      
      # if discrete, use the discrete initialization
      
      #cpt <- new("CondProbTable",
      #           factors = elim.tree[[this.tree]]@members,
      #           table = model$cptable,
      #           logprob = model$probs
      # )
      
      # elim.tree[[this.tree]]@cpt <- cpt
      
      
    } else {
      # if continuous, use the continuous initialization
      elim.tree <- InitCont(elim.tree, this.tree, dag, model$bags)
    }
  }
  
  result <- list(nodes=graph::nodes(dag), 
                 node.class=node.class,
                 dag=dag, 
                 tree.graph=tree, 
                 cluster.tree=elim.tree,
                 assignment=asgn,
                 propagated=FALSE,
                 absorbed.variables=character(0),
                 absorbed.values=list(),
                 absorbed.soft.variables=character(0),
                 absorbed.soft.values=list())
  class(result) <- "clustertree"
  
  return(result)
}




