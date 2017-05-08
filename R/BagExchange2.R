

BagExchange.2 <- function (bag1, bag2) {
  bag.post <- bag1
  bag.lp <- bag2
  
  conf.1 <- cond.extract(bag.post)
  conf.2 <- cond.extract(bag.lp)
  
  condv.1 <- bag.post[[1]]@conditionvals
  condv.2 <- bag.lp[[1]]@conditionvals 
  
  com.v <- intersect(condv.1, condv.2)
  
  if(length(com.v)==0) {
    
    exchg <- CompleteExchange(bag.post, bag.lp)
    output.lp <- exchg[[1]]
    output.post <- exchg[[2]]
    output <- list()
    output[[1]] <- output.lp
    output[[2]] <- output.post
    return(output)
  }
  
  
  confsub.1 <- conf.1[com.v]
  confsub.2 <- conf.2[com.v]
  
  conf.char.1 <- apply(confsub.1, 1, paste0, collapse="%")
  conf.char.2 <- apply(confsub.2, 1, paste0, collapse="%")
  
  #conf.char.1.s <- sort(conf.char.1)
  #conf.char.2.s <- sort(conf.char.2)
  
  conf.uniq <- unique(conf.char.1)
  
  smbags.post <- list()
  smbags.lp <- list()
  
  output.post <- list()
  output.lp <- list()
  
  for (i in 1:length(conf.uniq)) {
    bg.1 <- which(conf.char.1 == conf.uniq[i])
    smbags.post[[i]] <- bag.post[bg.1]
    bg.2 <- which(conf.char.2 == conf.uniq[i])
    smbags.lp[[i]] <- bag.lp[bg.2]
    
    exchg <- CompleteExchange(smbags.post[[i]], smbags.lp[[i]])
    output.lp <- c(output.lp, exchg[[1]])
    output.post <- c(output.post, exchg[[2]])
  }
  
  output <- list()
  output[[1]] <- output.lp
  output[[2]] <- output.post
    
  return(output)  
}


###########################################
## Complete exhange two bags of potentials without checking conditions
###########################################

CompleteExchange <- function(bag1, bag2){
  output <- newbag1 <- newbag2 <- list()
  for (i in 1:length(bag1)) {
    for (j in 1:length(bag2)) {
      lp1 <- bag1[[i]]
      lp2 <- bag2[[j]]
      #if (is.Exchangeable(lp1, lp2)) {
      newlps <- Exchange(lp1, lp2)
      newbag1 <- c(newbag1, newlps[[1]])
      newbag2 <- c(newbag2, newlps[[2]])
      #}
    }
  }
  output[[1]] <- newbag1
  output[[2]] <- newbag2
  return(output)
}

######## helper functions ################

extract.single <- function(lp){return(lp@conditionvalues)}


cond.extract <- function(bag) {
  # bag <- bag1
  
  if(length(bag)>0) {
    cond.vars <- bag[[1]]@conditionvals
  } else {
    return()
  }
  
  #conds <- vector("list", length=length(bag))
  #for (i in 1:length(bag)){
  #  conds[[i]] <- bag[[i]]@conditionvalues
  #}
  
  conds <- lapply(bag, extract.single)
  mat <- matrix(unlist(conds), nrow=length(bag), byrow=TRUE)
  
  conds.df <- data.frame(mat, stringsAsFactors=FALSE)
  colnames(conds.df) <- cond.vars
  return(conds.df)
}

