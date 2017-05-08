
PlotPosteriorDiscrete <- function(posteriors, group=NULL) {
  
  if(is.null(group)){
    nms <- names(posteriors)
  } else {
    nms <- group
  }
  
  ncl <- length(posteriors[[1]])
  if (ncl<3) ncl <- 3
  
  barplot(t(t(posteriors[[1]])), width=0.5, space=0.2, main=nms,
          col=brewer.pal(ncl, "Blues"), xlab=names(posteriors),
          beside=TRUE, names.arg=names(posteriors[[1]]))
  
}


PlotPosteriorDiscrete.old <- function(posteriors) {
  
  n <- length(posteriors)
  nms <- names(posteriors)
  
  if(n==1){
    barplot(t(t(posteriors[[1]])), width=0.5, space=0.2, main=nms,
            col=brewer.pal(length(posteriors[[1]]), "Blues"), 
            beside=TRUE, names.arg=names(posteriors[[1]]))
  }
  
  if (n>1){
    barplot(t(t(posteriors[[1]])), xlim=c(0,(n-1)*0.75+0.5), width=0.5, space=0, 
            col=brewer.pal(length(posteriors[[1]]), "Blues"), names.arg=nms[1])
    for (i in 2:n){
      barplot(t(t(posteriors[[i]])), add=TRUE, width=0.5, 
              offset=0, space=1.5*i-1.5, axes=FALSE,  col=brewer.pal(length(posteriors[[i]]), "Blues"),
              names.arg=nms[i])
    }
  }
  
}