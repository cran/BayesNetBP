#' Obtain marginal distributions
#'
#' Get the marginal distributions of multiple variables
#'
#' @details Get the marginal distributions of multiple variables. The function \code{Marginals}
#' returns a \code{list} of marginal distributions. The marginal distribution of a discrete variable 
#' is a named vector of probabilities. Meanwhile, the marginal distributions of 
#' continous variables in a CG-BN model are mixtures of Gaussian distributions. 
#' To fully represent this information, the marginal of a continuous variable is represented by 
#' a \code{data.frame} with three columns to specify 
#' parameters for each Gaussian distribution in the mixture, which are
#' 
#' \describe{
#'  \item{\code{mean}}{the mean value of a Gaussian distribution.}
#'  \item{\code{sd}}{the standard deviation of a Gaussian distribution.}
#'  \item{\code{n}}{the number of Gaussian mixtures}
#' }
#'
#' @param tree a \code{clustertree} object
#' @param vars a \code{vector} of variables for query of marginal distributions
#' 
#' @return 
#' 
#' \describe{
#'  \item{\code{marginals}}{a \code{list} of marginal distributions}
#'  \item{\code{types}}{a named \code{vector} indicating the types of the variables whose
#'  marginals are queried: \code{TRUE} for discrete, \code{FALSE} for continuous.}
#' }
#'
#' @author Han Yu
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. 
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
#' marg <- Marginals(tree.post, c("HDL", "Ppap2a"))
#' marg$marginals$HDL
#' head(marg$marginals$Ppap2a)
#' 
#' @seealso \code{\link{PlotMarginals}} for visualization of the marginal distributions,
#' \code{\link{SummaryMarginals}} for summarization of the marginal distributions of 
#' continuous variables.
#' 
#' @export

Marginals <- function(tree, vars) {
  
  object <- tree
  
  if(sum(vars %in% object$absorbed.variables)!=0) {
    var.in <- vars[vars %in% object$absorbed.variables]
    msg1 <- paste0(var.in, collapse=", ")
    stop(paste0(msg1, " is/are already observed."))
  }
  
  tree <- object$cluster.tree
  node.class <- object$node.class
  marginal.types <- node.class[vars]
  
  result <- list()
  for (i in 1:length(vars)) {
    var <- vars[i]
    if (node.class[[var]]) {
      result[[i]] <- DiscreteMarginal(tree, var)
    } else {
      result[[i]] <- PushMarginal(tree, var)
    }
  }
  names(result) <- vars
  
  output <- list()
  output$marginals <- result
  output$types <- marginal.types
  
  return(output)
}

####################

DiscreteMarginal <- function(post.tree, var){
  
  for(i in 1:length(post.tree)){
    if(var %in% post.tree[[i]]@members & post.tree[[i]]@discrete) {
      this.cluster <- post.tree[[i]]@name
      break
    }
  }
  
  df <- data.frame(post.tree[[this.cluster]]@joint$cpt, prob=post.tree[[this.cluster]]@joint$prob)
  fml <- as.formula(paste0("prob~", var))
  df.marg <- summaryBy(fml, df, FUN=sum)
  result <- df.marg[,2]
  names(result) <- df.marg[,1]
  
  return(result)
}

####################

PushMarginal <- function(post.tree, var){
  
  # Step 1
  n <- which(names(post.tree)==var)
  post.tree[[n]]@postbag <- post.tree[[n]]@lppotential
  
  
  # Step 3
  while((length(post.tree[[n]]@parent)!=0) & 
        !post.tree[[post.tree[[n]]@parent]]@discrete){
    post.tree[[post.tree[[n]]@parent]]@postbag <- post.tree[[n]]@postbag
    
    #########################################
    flag <- post.tree[[post.tree[[n]]@parent]]@activeflag
    if (length(post.tree[[post.tree[[n]]@parent]]@lppotential)==0){
      flag <- FALSE
    } else {
      lp.head <- post.tree[[post.tree[[n]]@parent]]@lppotential[[1]]@head
      postbag.tail <- post.tree[[post.tree[[n]]@parent]]@postbag[[1]]@tail
      if( !lp.head %in% postbag.tail){
        flag <- FALSE
      }
    }
    #########################################
    
    #if (post.tree[[post.tree[[n]]@parent]]@activeflag){
    if (flag) {
      newBag <- BagExchange.2(post.tree[[post.tree[[n]]@parent]]@postbag,
                              post.tree[[post.tree[[n]]@parent]]@lppotential)
      post.tree[[post.tree[[n]]@parent]]@postbag <- newBag[[2]]
      post.tree[[post.tree[[n]]@parent]]@lppotential <- newBag[[1]]
    }
    
    post.tree[[n]]@postbag <- list()
    
    n <- which(names(post.tree) == post.tree[[n]]@parent)
    
  }
  
  
  parent.cpt <- post.tree[[post.tree[[n]]@parent]]@joint$cpt #cpt@table
  post.parent.prob <- post.tree[[post.tree[[n]]@parent]]@joint$prob #cpt@logprob
  
  prob <- mu <- sd <- c()
  k <- 1
  
  post.vars <- post.tree[[n]]@postbag[[1]]@conditionvals
  cpt.sub <- parent.cpt[post.vars]
  
  cpt.char <- apply(cpt.sub, 1, paste0, collapse="%")
  
  
  for (i in 1:length(post.tree[[n]]@postbag)) {
    post.pot <- post.tree[[n]]@postbag[[i]]
    this.char <- paste0(post.pot@conditionvalues, collapse="%")
    
    pos <- which(cpt.char==this.char)
    this.prob <- sum(post.parent.prob[pos])
    
    prob[i] <- this.prob
    mu[i] <- post.pot@const
    sd[i] <- sqrt(post.pot@sigma)
  }
  result <- data.frame(prob, mu, sd)
  
  return(result)
}

