## Methods for crlmm

setGeneric("calls", function(x) standardGeneric("calls"))
setMethod("calls", "crlmmSet", function(x) assayData(x)$calls)

setGeneric("confs", function(x) standardGeneric("confs"))
setMethod("confs", "crlmmSet", function(x) assayData(x)$confs)

setGeneric("list2crlmmSet", function(x) standardGeneric("list2crlmmSet"))
setMethod("list2crlmmSet", "list",
          function(x){
            pd <- data.frame(SNR=x[["SNR"]],
                             row.names=colnames(x[["calls"]]))
            pdv <- data.frame(labelDescription=c("Signal-to-noise Ratio"),
                              row.names=c("SNR"))
            fd <- data.frame(SNPQC=x[["SNPQC"]],
                             row.names=rownames(x[["calls"]]))
            fdv <- data.frame(labelDescription=c("SNP Quality Score"),
                              row.names=c("SNPQC"))
            new("crlmmSet",
                assayData=assayDataNew("lockedEnvironment",
                  calls=x[["calls"]], confs=x[["confs"]]),
                phenoData=new("AnnotatedDataFrame",
                  data=pd, varMetadata=pdv),
                featureData=new("AnnotatedDataFrame",
                  data=fd, varMetadata=fdv))
          })
