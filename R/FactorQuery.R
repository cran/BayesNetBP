#' Queries of discrete variable distributions 
#'
#' Obtain the joint, marginal, and conditional distributions of discrete variables
#'
#' @details Query the joint distribution of any combination of discrete variables when 
#' mode is "joint", or conditional distribution of a discrete variable. The mode "list"
#' return a \code{list} of variable combinations, such that joint distributions of any subset 
#' of them are ready for extraction. Queries outside this list are also supported but may 
#' take longer computing time. This function will also return marginal distribution if only
#' one variable is queried.
#' 
#' @param tree a \code{clustertree} object
#' @param vars the variables to be queried
#' @param mode type of desired distribution
#' @return \code{data.frame} object specifying a joint or conditional distribution.
#' 
#' @author Han Yu
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. 
#' 
#' @importFrom igraph neighbors all_simple_paths induced_subgraph
#' 
#' @examples 
#' 
#' data(chest)
#' dag <- chest$dag
#' node.class <- rep(TRUE, length(dag@nodes))
#' names(node.class) <- dag@nodes
#' cst <- ClusterTreeCompile(dag, node.class)
#' models <- LocalModelCompile(chest$data, dag, node.class)
#' tree.init <- ElimTreeInitialize(tree=cst$tree.graph, 
#'                                 dag=cst$dag, 
#'                                 model=models, 
#'                                 node.sets=cst$cluster.sets, 
#'                                 node.class=cst$node.class)
#' tree.init.p <- PropagateDBN(tree.init)
#' # get joint distribution
#' FactorQuery(tree=tree.init.p, vars=c("tub", "xray", "dysp", "asia"), mode="joint")
#' 
#' # get joint distribution
#' FactorQuery(tree=tree.init.p, vars=c("xray"), mode="conditional")
#' 
#' @export

FactorQuery <- function(tree, vars=c(), mode=c("joint", "conditional", "list")) {
  
  if(sum(vars %in% tree$absorbed.variables)!=0) {
    var.in <- vars[vars %in% tree$absorbed.variables]
    msg1 <- paste0(var.in, collapse=", ")
    stop(paste0(msg1, " is/are already observed."))
  }
  
  cstree <- tree$cluster.tree
  dag <- tree$dag
  dag.graph <- igraph.from.graphNEL(dag)
  
  if (mode=="list") {
    result <- list()
    j <- 1
    for (i in 1:length(cstree)){
      if(cstree[[i]]@discrete & is.subset(vars, cstree[[i]]@members)) {
        result[[j]] <- cstree[[i]]@members
        j <- j+1
      }
    }
    return(result)
  }
  
  if (mode=="conditional") {
    if (length(vars)!=1) {
      stop("If mode is conditional, vars must be a single variable.")
    }
    parents <- names(neighbors(dag.graph, vars, mode="in"))
    parents <- setdiff(parents, tree$absorbed.variables)
    allvar <- c(vars, parents)
    for (i in 1:length(cstree)){
      if(cstree[[i]]@discrete & is.subset(allvar, cstree[[i]]@members)) {
        result <- cstree[[i]]@joint
        result <- marginalize.discrete(result, allvar)
        result <- conditional(result, parents)
        result.df <- data.frame(result$cpt, prob=result$prob)
        break
      }
    }
    rownames(result.df) <- NULL
    return(result.df)
  }
  
  if (mode=="joint") {
    
    in.cluster <- FALSE
    
    for (i in 1:length(cstree)){
      if(cstree[[i]]@discrete & is.subset(vars, cstree[[i]]@members)) {
        result <- cstree[[i]]@joint
        result <- marginalize.discrete(result, vars)
        result.df <- data.frame(result$cpt, prob=result$prob)
        in.cluster <- TRUE
        break
      }
    }
    
    if (!in.cluster){
      result <- query.ooc(tree, vars)
      result.df <- data.frame(result$cpt, prob=result$prob)
    }
    rownames(result.df) <- NULL
    return(result.df)
  }
}