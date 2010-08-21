##setGeneric("batch", function(object) standardGeneric("batch"))
##setGeneric("getParam", function(object, name, batch) standardGeneric("getParam"))
setGeneric("cnIndex", function(object) standardGeneric("cnIndex"))
setGeneric("cnNames", function(object) standardGeneric("cnNames"))
setGeneric("computeCopynumber", function(object, ...) standardGeneric("computeCopynumber"))
##setGeneric("pr", function(object, name, batch, value) standardGeneric("pr"))
setGeneric("snpIndex", function(object) standardGeneric("snpIndex"))
setGeneric("snpNames", function(object) standardGeneric("snpNames"))

setGeneric("totalCopyNumber", function(object,...) standardGeneric("totalCopyNumber"))

setGeneric("CA", function(object, i, j, ...) standardGeneric("CA"))
setGeneric("CB", function(object, i, j, ...) standardGeneric("CB"))
setGeneric("totalCopyNumber", function(object, i, j, ...) standardGeneric("totalCopyNumber"))


## The generics below are for internal use with copy number methods
## If we keep them in oligoClasses, we need to export and document
setGeneric("nuA", function(object) standardGeneric("nuA"))
setGeneric("nuB", function(object) standardGeneric("nuB"))
setGeneric("phiA", function(object) standardGeneric("phiA"))
setGeneric("phiB", function(object) standardGeneric("phiB"))
setGeneric("sigma2A", function(object) standardGeneric("sigma2A"))
setGeneric("sigma2B", function(object) standardGeneric("sigma2B"))
setGeneric("tau2A", function(object) standardGeneric("tau2A"))
setGeneric("tau2B", function(object) standardGeneric("tau2B"))
setGeneric("corrAA", function(object) standardGeneric("corrAA"))
setGeneric("corrBB", function(object) standardGeneric("corrBB"))
setGeneric("corrAB", function(object) standardGeneric("corrAB"))
setGeneric("nuA<-", function(object, value) standardGeneric("nuA<-"))
setGeneric("nuB<-", function(object, value) standardGeneric("nuB<-"))
setGeneric("phiA<-", function(object, value) standardGeneric("phiA<-"))
setGeneric("phiB<-", function(object, value) standardGeneric("phiB<-"))
setGeneric("sigma2A<-", function(object, value) standardGeneric("sigma2A<-"))
setGeneric("sigma2B<-", function(object, value) standardGeneric("sigma2B<-"))
setGeneric("tau2A<-", function(object, value) standardGeneric("tau2A<-"))
setGeneric("tau2B<-", function(object, value) standardGeneric("tau2B<-"))
setGeneric("corrAA<-", function(object, value) standardGeneric("corrAA<-"))
setGeneric("corrAB<-", function(object, value) standardGeneric("corrAB<-"))
setGeneric("corrBB<-", function(object, value) standardGeneric("corrBB<-"))


