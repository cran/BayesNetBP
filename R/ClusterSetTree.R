
#' An S4 class to represent a cluster tree node.
#'
#' 

setClass("ClusterSetTree",
         slots = list(name = "character", 
                      members = "vector",
                      discrete = "logical",
                      index = "numeric",
                      parent = "character",
                      children = "vector",
                      lppotential = "list",
                      postbag = "list",
                      cpt = "CondProbTable",
                      activeflag = "logical",
                      logweighttable = "list",
                      joint = "list"
         )
)








