setGeneric("batch", function(object) standardGeneric("batch"))
setGeneric("getParam", function(object, name, batch) standardGeneric("getParam"))
setGeneric("cnIndex", function(object) standardGeneric("cnIndex"))
setGeneric("cnNames", function(object) standardGeneric("cnNames"))
setGeneric("computeCopynumber", function(object, ...) standardGeneric("computeCopynumber"))
setGeneric("pr", function(object, name, batch, value) standardGeneric("pr"))
setGeneric("snpIndex", function(object) standardGeneric("snpIndex"))
setGeneric("snpNames", function(object) standardGeneric("snpNames"))
setGeneric("lM", function(object) standardGeneric("lM"))
setGeneric("lM<-", function(object, value) standardGeneric("lM<-"))
setGeneric("totalCopyNumber", function(object,...) standardGeneric("totalCopyNumber"))

setGeneric("corr", function(object, allele) standardGeneric("corr"))
setGeneric("nu", function(object, allele) standardGeneric("nu"))
setGeneric("phi", function(object, allele) standardGeneric("phi"))
setGeneric("sigma2", function(object, allele) standardGeneric("sigma2"))
setGeneric("tau2", function(object, allele) standardGeneric("tau2"))

setGeneric("CA", function(object, i, j, ...) standardGeneric("CA"))
setGeneric("CB", function(object, i, j, ...) standardGeneric("CB"))
setGeneric("totalCopyNumber", function(object, i, j, ...) standardGeneric("totalCopyNumber"))


