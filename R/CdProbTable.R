

setClass("CondProbTable",
         slots = list(node = "character", 
                      factors = "vector",
                      table = "data.frame",
                      logprob = "vector",
                      values = "list",
                      probs = "list"
         )
)








