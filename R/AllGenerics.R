##setGeneric("A<-", function(object, value) standardGeneric("A<-"))
##setGeneric("B<-", function(object, value) standardGeneric("B<-"))

setGeneric("getParam", function(object, name, ...) standardGeneric("getParam"))
setGeneric("cnIndex", function(object) standardGeneric("cnIndex"))
setGeneric("cnNames", function(object) standardGeneric("cnNames"))
setGeneric("computeCopynumber", function(object, cnOptions) standardGeneric("computeCopynumber"))
setGeneric("pr", function(object, name, batch, value) standardGeneric("pr"))
setGeneric("snpIndex", function(object) standardGeneric("snpIndex"))
setGeneric("snpNames", function(object) standardGeneric("snpNames"))
##setGeneric("splitByChromosome", function(object, ...) standardGeneric("splitByChromosome"))

