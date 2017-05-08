#' A Cluster Tree class
#'
#' The \code{clustertree} object is the computational object for belief propagation.
#' 
#' @name clustertree
#' @details A \code{clustertree} object can be obtained by the  \code{ElimTreeInitialize}, 
#' and it contains necessary information for belief propagation, such as
#' the topological inforamtion of Bayesian network and cluster tree, as well as the local models. This
#' object comprises the following elements,
#' \describe{
#'   \item{\code{nodes}}{A \code{list} of the node or variable names.}
#'   \item{\code{node.class}}{A named \code{vector} of logical values indicating whether each node is continuous or discrete. }
#'   \item{\code{dag}}{A \code{graphNEL} object of the Bayesian network structure.}
#'   \item{\code{tree.graph}}{A \code{graphNEL} object of the semi-elimination tree structure. }
#'   \item{\code{cluster.tree}}{ A \code{list} storing the cluster members and corresponding local models as linear predictor
#'    potentials, conditional probability tables, and joint distribution tables.}
#'   \item{\code{assignment}}{ A \code{list} indicating the cluster assignment of the discrete factors. }
#'   \item{\code{propagated}}{ A \code{logical} value indicating whether the discrete compartment has been propagated. }
#'   \item{\code{absorbed.variables}}{ A \code{vector} of characters indicating variables observed with hard evidence. }
#'   \item{\code{absorbed.values}}{ A \code{list} indicating the values of the variables observed with hard evidence. }
#'   \item{\code{absorbed.soft.variables}}{ A \code{vector} of characters indicating variables observed with soft or 
#'   likelihood evidence. }
#'   \item{\code{absorbed.soft.values}}{ A \code{list} of the likelihoods of the soft or likelihood evidence. }
#' } 
#' 
#' @seealso \code{\link{ElimTreeInitialize}}
#'

NULL