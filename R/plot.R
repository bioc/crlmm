plotSamples = function(cnSet, col=1, offset=0, xlim=c(9,16), ylim=c(-5,5), verbose=FALSE, sample=100000, seed=1, type="smoothScatter"){
    if(missing(col)){
        stop("col is missing, please specify which sample/(s) to plot")
    }
    set.seed(seed)
    row = sample(nrow(cnSet), sample)
    if(verbose)  
      message("compute M values")        
    M = computeLogRatio(cnSet, offset=offset, verbose=verbose, row=row, col=col)
    if(type=="beanplot") {
      sel = rowSums(is.finite(M[,]) & !is.na(M[,]))==ncol(M)
      beanplot(as.vector(M[sel,])~factor(rep(col, each=sum(sel))), ylab="M", xlab="Sample", ylim=ylim,
               beanlines="median", what=c(0,1,1,0), col="grey", border="black") #c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", 
      abline(h=0, col="gray", lty=2)
      delete(M); rm(M)
    }
    if(type=="smoothScatter") {
      if(verbose)
        message("compute S values")
      S = computeAverageLogIntensity(cnSet, offset=offset, verbose=verbose, row=row, col=col)
      for(i in 1:length(col)) {
        smoothScatter(S[,i], M[,i], main=colnames(M)[i], ylab="M", xlab="S", xlim=xlim, ylim=ylim)
        abline(h=0, col="gray", lty=2)
      }
    delete(S, M); rm(S, M)
    }
}


plotSNPs = function(cnSet, row=1, offset=0, xlim=c(9,16), ylim=c(-5,5), verbose=FALSE){
    if(missing(row)){
        stop("row is missing, please specify which SNP/(s) to plot")
    }
    if(verbose)  
      message("compute M values")        
    M = computeLogRatio(cnSet, offset=offset, verbose=verbose, row=row)
    if(verbose)
      message("compute S values")   
    S = computeAverageLogIntensity(cnSet, offset=offset, verbose=verbose, row=row)
    if(verbose)
      message("get genotype calls")
    for(i in 1:length(row)){
       plot(S[i,], M[i,], xlab="S", ylab="M", xlim=xlim, ylim=ylim,
            main=rownames(M)[i], pch=19, col=calls(cnSet)[row[i],])
       abline(h=0, col="gray", lty=2)
    }
    rm(S, M)
}
