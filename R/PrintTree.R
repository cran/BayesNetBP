#' Print the cluster tree
#'
#' Print a \code{clustertree} object in a readable format
#' 
#' @details Print a \code{clustertree} object in a readable format by outputing cluster members
#' and linear potentials assigned to each cluster.
#'
#' @param tree a \code{clustertree} object
#' 
#' @author Han Yu
#' 
#' @examples 
#' data(toytree)
#' PrintTree(toytree)
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550.  
#' 
#' @seealso \code{\link{clustertree}}
#' 
#' @export

PrintTree <- function(tree){
  tree.print <- tree$cluster.tree
  for (j in 1:length(tree.print)) {
    cat(tree.print[[j]]@name, ":", tree.print[[j]]@members, "\n")
    cat("LPPotential: \n")
    
    if (length(tree.print[[j]]@lppotential)>0){
      for (k in 1:length(tree.print[[j]]@lppotential)){
        cat(tree.print[[j]]@lppotential[[k]]@conditionvals, "\t")
        cat(tree.print[[j]]@lppotential[[k]]@conditionvalues, "\t")
        cat(tree.print[[j]]@lppotential[[k]]@head, "|")
        cat(tree.print[[j]]@lppotential[[k]]@tail, ": ")
        cat(tree.print[[j]]@lppotential[[k]]@params, 
            "; constant=", tree.print[[j]]@lppotential[[k]]@const,
            "; sigma=", tree.print[[j]]@lppotential[[k]]@sigma, "\n")
      }
    }
    
    cat("Postbag: \n")
    if (length(tree.print[[j]]@postbag)>0){
      for (k in 1:length(tree.print[[j]]@postbag)){
        cat(tree.print[[j]]@postbag[[k]]@conditionvals, "\t")
        cat(tree.print[[j]]@postbag[[k]]@conditionvalues, "\t")
        cat(tree.print[[j]]@postbag[[k]]@head, "|")
        cat(tree.print[[j]]@postbag[[k]]@tail, ": ")
        cat(tree.print[[j]]@postbag[[k]]@params, 
            "; constant=", tree.print[[j]]@postbag[[k]]@const,
            "; sigma=", tree.print[[j]]@postbag[[k]]@sigma, "\n")
      }
    }
    cat("__________________________________\n \n")
  }
}