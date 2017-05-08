
ModelCompileData <- function(data, dag, node.class) {
  
  # qtl <- small_liver; qtl.fit <- out.liver; pheno.cols=2:11
  
  ###########
  
  # dag, data.frame, node.class
  # data <- dat
  nodes <- dag@nodes
  dag.graph <- igraph.from.graphNEL(dag)
  
  value.list <- list()
  discrete.nodes <- nodes[node.class]
  continuous.nodes <- nodes[!node.class]
  
  df <- data
  
  dat.complete.0 <- df[complete.cases(df),]
  
  ######################
  ## discrete part starts
  ######################
  cpt.pots <- list()
  
  if (length(discrete.nodes)>0){
    
    for (i in 1:length(discrete.nodes)) {
      value.list[[i]] <- sort(unique(dat.complete.0[[discrete.nodes[i]]]))
    }
    names(value.list) <- discrete.nodes
    
    # source("gRain_test_functions.R")
    
    
    
    for (i in 1:length(discrete.nodes)) {
      # i <- 2
      this.node <- discrete.nodes[i]
      this.parents <- names(neighbors(dag.graph, this.node, "in"))
      all.nodes <- c(this.node, this.parents)
      n.nodes <- length(all.nodes)
      this.cpt <- expand.grid(value.list[all.nodes])
      this.df <- df[all.nodes]
      
      tab <- as.data.frame(xtabs(~., this.df))
      pot.joint <- list(cpt=tab[1:n.nodes], prob=tab$Freq/sum(tab$Freq))
      
      if(length(this.parents)==0) {
        pot <- pot.joint
      } else {
        pot <- conditional(pot.joint, this.parents)
      }
      
      cpt.pots[[i]] <- pot
    }
    
    names(cpt.pots) <- discrete.nodes
    
  }
  
  ######################
  ## discrete part ends
  ######################
  
  
  ######################
  ## continuous part starts
  ######################
  
  bags <- list()
  
  if (length(continuous.nodes)>0){
    
    bags <- vector("list", length(continuous.nodes))
    
    for (i in 1:length(continuous.nodes)){
      # i <- 1
      
      this.bag <- vector("list", 1000)
      k <- 1
      
      this.node <- continuous.nodes[i]
      this.parents <- names(neighbors(dag.graph, this.node, "in"))
      
      dat.complete <- df[c(this.node, this.parents)]
      dat.complete <- dat.complete[complete.cases(dat.complete),]
      
      this.classes <- node.class[this.parents]
      
      ######################
      
      discrete.parents <- this.parents[which(this.classes)]
      continuous.parents <- this.parents[which(!this.classes)]
      
      ######################
      
      if(length(discrete.parents)==0 & length(continuous.parents)==0){
        y <- dat.complete[[this.node]]
        lp <- new("LPPotential", 
                  head = this.node,
                  const = mean(y),
                  sigma = var(y)
        )
        this.bag[[k]] <- lp
        k <- k+1
        this.bag <- Filter(Negate(is.null), this.bag)
        bags[[i]] <- this.bag
        next
      }
      
      ######################
      
      if(length(discrete.parents)==0){
        
        df.2 <- dat.complete[c(this.node, continuous.parents)]
        df.sub <- df.2
        form.str <- paste0(this.node, "~.")
        form <- as.formula(form.str)
        lm.fit <- lm(form, df.sub)
        coefs <- coef(lm.fit)
        lp <- new("LPPotential", 
                  head = this.node, 
                  tail = continuous.parents,
                  params = coefs[2:length(coefs)],
                  const = coefs[1],
                  sigma = summary(lm.fit)$sigma^2)
        this.bag[[k]] <- lp
        k <- k+1
        this.bag <- Filter(Negate(is.null), this.bag)
        bags[[i]] <- this.bag
        next
      }
      
      ######################
      
      this.disc.vals <- value.list[discrete.parents]
      this.all.combs <- expand.grid(this.disc.vals, stringsAsFactors=FALSE)
      comb.val.list <- apply(this.all.combs, 1, paste0, collapse="%")
      df.1 <- dat.complete[discrete.parents]  ## mark
      df.comb <- apply(df.1, 1, paste0, collapse="%")
      
      
      if(length(continuous.parents)==0){
        for (j in 1:length(comb.val.list)) {
          # j <- 1
          sub.ind <- which(df.comb==comb.val.list[j])
          df.2 <- dat.complete[[this.node]]
          y <- df.2[sub.ind]
          lp <- new("LPPotential", 
                    head = this.node, 
                    const = mean(y),
                    sigma = var(y),
                    conditionvals = discrete.parents,
                    conditionvalues = as.vector(this.all.combs[j,], 
                                                mode="character") )
          this.bag[[k]] <- lp
          k <- k+1
        }
        this.bag <- Filter(Negate(is.null), this.bag)
        bags[[i]] <- this.bag
        next
      }
      
      
      
      for (j in 1:length(comb.val.list)) {
        # j <- 1
        sub.ind <- which(df.comb==comb.val.list[j])
        df.2 <- dat.complete[c(this.node, continuous.parents)]
        df.sub <- df.2[sub.ind,]
        form.str <- paste0(this.node, "~.")
        form <- as.formula(form.str)
        lm.fit <- lm(form, df.sub)
        coefs <- coef(lm.fit)
        lp <- new("LPPotential", 
                  head = this.node, 
                  tail = continuous.parents,
                  params = coefs[2:length(coefs)],
                  const = coefs[1],
                  sigma = summary(lm.fit)$sigma^2,
                  conditionvals = discrete.parents,
                  conditionvalues = as.vector(this.all.combs[j,], 
                                              mode="character") )
        this.bag[[k]] <- lp
        k <- k+1
      }
      this.bag <- Filter(Negate(is.null), this.bag)
      bags[[i]] <- this.bag
    }
    
    names(bags) <- continuous.nodes
    
  }
  
  ######################
  ## continuous part ends
  ######################
  
  result <- list(pots = cpt.pots,
                 bags = bags)
  
  
  return(result)
  
}

##################################################
## helper function
##################################################

loci.loc <- function(locus) {
  split <- strsplit(locus, split='_')
  chr <- c()
  location <- c()
  for (i in 1:length(split)) {
    chr[i] <- substring(split[[i]][1], 4)
    location[i] <- split[[i]][2]
  }
  
  result <- data.frame(chr, location)
  result$chr <- as.character(result$chr)
  result$location <- as.numeric(as.character(result$location))
  
  return(result)
}


##########################################
## Initialize 
##########################################

init.cont <- function(elim.tree, node, dag, lp.bags){
  this.node <- node
  cs.tree <- elim.tree
  parent.nodes <- parents(this.node, dag)
  
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