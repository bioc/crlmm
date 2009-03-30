# function below works OK provided all .idat files are in the current working directory
# - could add an option to allow files in Illumina directory structure to be handled
# or to use the optional 'Path' column in sampleSheet
# - there is a lot of header information that is currently discarded - could try and store this somewhere in the resulting NChannelSet

readIdatFiles = function(sampleSheet=NULL, arrayNames=NULL, ids=NULL, path=".",
                                arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
                                highDensity=FALSE, sep="_", fileExt=list(green="Grn.idat", red="Red.idat")) {
       if(!is.null(sampleSheet)) { # get array info from Illumina's sample sheet
         arrayNames=NULL
         if(!is.null(arrayInfoColNames$barcode) && (arrayInfoColNames$barcode %in% colnames(sampleSheet))) {
           barcode = sampleSheet[,arrayInfoColNames$barcode]
           arrayNames=barcode
         }
         if(!is.null(arrayInfoColNames$position) && (arrayInfoColNames$position %in% colnames(sampleSheet))) {  
           position = sampleSheet[,arrayInfoColNames$position]
           if(is.null(arrayNames))
             arrayNames=position
           else
             arrayNames = paste(arrayNames, position, sep=sep)
           if(highDensity) {
             hdExt = list(A="R01C01", B="R01C02", C="R02C01", D="R02C02")
             for(i in names(hdExt))
               arrayNames = sub(paste(sep, i, sep=""), paste(sep, hdExt[[i]], sep=""), arrayNames)
            }
         }
#         pd = new("AnnotatedDataFrame", data = sampleSheet)
       }
       if(is.null(arrayNames)) {
         arrayNames = gsub(paste(sep, fileExt$green, sep=""), "", dir(pattern=fileExt$green, path=path))
         if(!is.null(sampleSheet)) {
           sampleSheet=NULL
           cat("Could not find required info in \'sampleSheet\' - ignoring.  Check \'sampleSheet\' and/or \'arrayInfoColNames\'\n")
         }
       }
       grnfiles = paste(arrayNames, fileExt$green, sep=sep)
       redfiles = paste(arrayNames, fileExt$red, sep=sep)
       if(length(grnfiles)==0 || length(redfiles)==0)
         stop("Cannot find .idat files")
       if(length(grnfiles)!=length(redfiles))
         stop("Cannot find matching .idat files")
       if(!all(c(redfiles,grnfiles) %in% dir(path=path)))
         stop("Missing .idat files: red\n", paste(redfiles[!(redfiles %in% dir(path=path))], sep=" "), "\n green\n",
                 paste(grnfiles[!(grnfiles %in% dir(path=path))], sep=" "))
       grnidats = file.path(path, grnfiles)
       redidats = file.path(path, redfiles)
       
       headerInfo = list(nProbes = rep(NA, length(arrayNames)),
                         Barcode = rep(NA, length(arrayNames)),
                         ChipType = rep(NA, length(arrayNames)),
                         Manifest = rep(NA, length(arrayNames)), # not sure about this one - sometimes blank
                         Position = rep(NA, length(arrayNames))) # this may also vary a bit
       # read in the data
       for(i in seq(along=arrayNames)) {
         cat("reading", arrayNames[i], "\t")
         idsG = idsR = G = R = NULL
#         if(grnfiles[i] %in% dir(path=path)) {
         cat(paste(sep, fileExt$green, sep=""), "\t")
         G = readIDAT(grnidats[i])
         idsG = rownames(G$Quants)
#         }
#         else
#           stop("Could not find ", grnidats[i])
         
#         if(redfiles[i] %in% dir(path=path)) {
          cat(paste(sep, fileExt$red, sep=""), "\n")
          R = readIDAT(redidats[i])
          idsR = rownames(R$Quants)
#         }
#         else
#           stop("Could not find ", redidats[i])
         
         headerInfo$nProbes[i] = G$nSNPsRead
         headerInfo$Barcode[i] = G$Barcode
         headerInfo$ChipType[i] = G$ChipType
         headerInfo$Manifest[i] = G$Unknown$MostlyNull
         headerInfo$Position[i] = G$Unknowns$MostlyA
         
         if(headerInfo$ChipType[i]!=headerInfo$ChipType[1] || headerInfo$Manifest[i]!=headerInfo$Manifest[1])
                  # || headerInfo$nProbes[i]!=headerInfo$nProbes[1] ## removed this condition as some arrays used the same manifest
                  # but differed by a few SNPs for some reason - most of the chip was the same though
           stop("Chips are not of all of the same type - please check your data")
         if(i==1) {
           if(is.null(ids) && !is.null(G))
             ids = idsG
           else
             stop("Could not find probe IDs")
           nprobes = length(ids)
           narrays = length(arrayNames)
           Gintens = Rintens = Gnobeads = Rnobeads = Gstderr = Rstderr = matrix(NA, nprobes, narrays)
           rownames(Gintens) = rownames(Rintens) = rownames(Gnobeads) = rownames(Rnobeads) = rownames(Gstderr) = rownames(Rstderr) = ids
           if(!is.null(sampleSheet))
             colnames(Gintens) = colnames(Rintens) =  colnames(Gnobeads) = colnames(Rnobeads) = colnames(Gstderr) = colnames(Rstderr) = sampleSheet$Sample_ID
           else
             colnames(Gintens) = colnames(Rintens) =  colnames(Gnobeads) = colnames(Rnobeads) = colnames(Gstderr) = colnames(Rstderr) = arrayNames
         }
         
         if(length(ids)==length(idsG)) {
           if(sum(ids==idsG)==nprobes) {
             Gintens[,i] = G$Quants[, "Mean"]
             Gnobeads[,i] = G$Quants[, "NBeads"]
             Gstderr[,i] = G$Quants[, "SD"]
           }
         }
         else {
           indG = match(ids, idsG)
           Gintens[,i] = G$Quants[indG, "Mean"]
           Gnobeads[,i] = G$Quants[indG, "NBeads"]
           Gstderr[,i] = G$Quants[indG, "SD"]
         }
         if(length(ids)==length(idsG)) {   
           if(sum(ids==idsR)==nprobes) {
             Rintens[,i] = R$Quants[ ,"Mean"]
             Rnobeads[,i] = R$Quants[ ,"NBeads"]
             Rstderr[,i] = R$Quants[ ,"SD"]
           }
         }
         else {
           indR = match(ids, idsR)
           Rintens[,i] = R$Quants[indR, "Mean"]
           Rnobeads[,i] = R$Quants[indR, "NBeads"]
           Rstderr[,i] = R$Quants[indR, "SD"]
         }
       }
       if(is.null(sampleSheet))
         pd = new("AnnotatedDataFrame", data = data.frame(Sample_ID=arrayNames))
       else
         pd = new("AnnotatedDataFrame", data = sampleSheet) 
       RG = new("NChannelSet",
                R=Rintens, G=Gintens, Rnb=Rnobeads, Gnb=Gnobeads,
                Rse=Rstderr, Gse=Gstderr, annotation=headerInfo$Manifest[1],
                phenoData = pd)
}


# the readIDAT() and readBPM() functions below were provided by Keith Baggerly, 27/8/2008
readIDAT <- function(idatFile){

  fileSize <- file.info(idatFile)$size

  tempCon <- file(idatFile,"rb")
  prefixCheck <- readChar(tempCon,4)
  if(prefixCheck != "IDAT"){

  }

  versionNumber <- readBin(tempCon, "integer", n=1, size=8, 
                           endian="little", signed=FALSE)
  
  nFields <- readBin(tempCon, "integer", n=1, size=4, 
                     endian="little", signed=FALSE)

  fields <- matrix(0,nFields,3);
  colnames(fields) <- c("Field Code", "Byte Offset", "Bytes")
  for(i1 in 1:nFields){
    fields[i1,"Field Code"] <- 
      readBin(tempCon, "integer", n=1, size=2, endian="little", signed=FALSE)
    fields[i1,"Byte Offset"] <- 
      readBin(tempCon, "integer", n=1, size=8, endian="little", signed=FALSE)
  }

  knownCodes <- c(1000, 102, 103, 104, 107, 200, 300, 400,
                  401, 402, 403, 404, 405, 406, 407, 408, 409)
  names(knownCodes) <- 
    c("nSNPsRead",  # 1000
      "IlluminaID", #  102
      "SD",         #  103
      "Mean",       #  104
      "NBeads",     #  107
      "MidBlock",   #  200
      "RunInfo",    #  300
      "RedGreen",   #  400
      "MostlyNull", #  401
      "Barcode",    #  402
      "ChipType",   #  403
      "MostlyA",    #  404
      "Unknown.1",  #  405
      "Unknown.2",  #  406
      "Unknown.3",  #  407
      "Unknown.4",  #  408
      "Unknown.5"   #  409
      )

  nNewFields <- 1
  rownames(fields) <- paste("Null", 1:nFields)
  for(i1 in 1:nFields){
    temp <- match(fields[i1,"Field Code"], knownCodes)
    if(!is.na(temp)){
      rownames(fields)[i1] <- names(knownCodes)[temp]
    }else{
      rownames(fields)[i1] <- paste("newField", nNewFields, sep=".")
      nNewFields <- nNewFields + 1
    }
  }

  seek(tempCon, fields["nSNPsRead", "Byte Offset"])
  nSNPsRead <- readBin(tempCon, "integer", n=1, size=4, 
                       endian="little", signed=FALSE)

  seek(tempCon, fields["IlluminaID", "Byte Offset"])
  IlluminaID <- readBin(tempCon, "integer", n=nSNPsRead, size=4, 
                       endian="little", signed=FALSE)

  seek(tempCon, fields["SD", "Byte Offset"])
  SD <- readBin(tempCon, "integer", n=nSNPsRead, size=2, 
                endian="little", signed=FALSE)

  seek(tempCon, fields["Mean", "Byte Offset"])
  Mean <- readBin(tempCon, "integer", n=nSNPsRead, size=2, 
                  endian="little", signed=FALSE)

  seek(tempCon, fields["NBeads", "Byte Offset"])
  NBeads <- readBin(tempCon, "integer", n=nSNPsRead, size=1, signed=FALSE)

  seek(tempCon, fields["MidBlock", "Byte Offset"])
  nMidBlockEntries <- readBin(tempCon, "integer", n=1, size=4, 
                              endian="little", signed=FALSE)
  MidBlock <- readBin(tempCon, "integer", n=nMidBlockEntries, size=4, 
                      endian="little", signed=FALSE)

  seek(tempCon, fields["RunInfo", "Byte Offset"])
  nRunInfoBlocks <- readBin(tempCon, "integer", n=1, size=4, 
                            endian="little", signed=FALSE)
  RunInfo <- matrix(NA, nRunInfoBlocks, 5)
  colnames(RunInfo) <- c("RunTime", "BlockType", "BlockPars", 
                         "BlockCode", "CodeVersion")
  for(i1 in 1:nRunInfoBlocks){
    for(i2 in 1:5){
      nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
      RunInfo[i1,i2] <- readChar(tempCon, nChars)
    }
  }

  seek(tempCon, fields["RedGreen", "Byte Offset"])
  RedGreen <- readBin(tempCon, "numeric", n=1, size=4, 
                      endian="little", signed=FALSE)
  #RedGreen <- readBin(tempCon, "integer", n=4, size=1, 
  #                    endian="little", signed=FALSE)

  seek(tempCon, fields["MostlyNull", "Byte Offset"])
  nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  MostlyNull <- readChar(tempCon, nChars) 

  seek(tempCon, fields["Barcode", "Byte Offset"])
  nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  Barcode <- readChar(tempCon, nChars) 

  seek(tempCon, fields["ChipType", "Byte Offset"])
  nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  ChipType <- readChar(tempCon, nChars) 

  seek(tempCon, fields["MostlyA", "Byte Offset"])
  nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  MostlyA <- readChar(tempCon, nChars) 

  seek(tempCon, fields["Unknown.1", "Byte Offset"])
  nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  Unknown.1 <- readChar(tempCon, nChars) 

  seek(tempCon, fields["Unknown.2", "Byte Offset"])
  nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  Unknown.2 <- readChar(tempCon, nChars) 

  seek(tempCon, fields["Unknown.3", "Byte Offset"])
  nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  Unknown.3 <- readChar(tempCon, nChars) 

  seek(tempCon, fields["Unknown.4", "Byte Offset"])
  nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  Unknown.4 <- readChar(tempCon, nChars) 

  seek(tempCon, fields["Unknown.5", "Byte Offset"])
  nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  Unknown.5 <- readChar(tempCon, nChars) 

  close(tempCon)

  Unknowns <- 
    list(MostlyNull=MostlyNull,
         MostlyA=MostlyA,
         Unknown.1=Unknown.1,
         Unknown.2=Unknown.2,
         Unknown.3=Unknown.3,
         Unknown.4=Unknown.4,
         Unknown.5=Unknown.5)

  Quants <- cbind(Mean, SD, NBeads)
  colnames(Quants) <- c("Mean", "SD", "NBeads")
  rownames(Quants) <- as.character(IlluminaID)

  idatValues <- 
    list(fileSize=fileSize, 
         versionNumber=versionNumber, 
         nFields=nFields, 
         fields=fields,
         nSNPsRead=nSNPsRead, 
         #IlluminaID=IlluminaID, 
         #SD=SD, 
         #Mean=Mean, 
         #NBeads=NBeads,
         Quants=Quants,
         MidBlock=MidBlock, 
         RunInfo=RunInfo, 
         RedGreen=RedGreen, 
         Barcode=Barcode, 
         ChipType=ChipType,
         Unknowns=Unknowns)

  idatValues

}

readBPM <- function(bpmFile){

  ## Reads and parses Illumina BPM files
  
  fileSize <- file.info(bpmFile)$size

  tempCon <- file(bpmFile,"rb")

  # The first few bytes of the egtFile are some type of
  # header, but there's no related byte offset information. 

  prefixCheck <- readChar(tempCon,3) ## should be "BPM"

  null.1 <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  ## should be 1  
  
  versionNumber <- 
    readBin(tempCon, "integer", n=1, size=4, endian="little", signed=FALSE)
  ## should be 4

  nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
  chipType <- readChar(tempCon, nChars)

  null.2 <- readBin(tempCon, "integer", n=2, size=1, signed=FALSE)

  csvLines <- readLines(tempCon, 22)

  entriesByteOffset <- seek(tempCon);
  nEntries <- readBin(tempCon, "integer", n=1, size=4,
                      endian="little", signed=FALSE)
    
  if(FALSE){
    
    snpIndexByteOffset <- seek(tempCon)
    snpIndex <- readBin(tempCon, "integer", n=nEntries, size=4,
                        endian="little", signed=FALSE)
    ## for the 1M array, these are simply in order from 1 to 1072820.
    
    snpNamesByteOffset <- seek(tempCon)
    snpNames <- rep("A", nEntries)
    for(i1 in 1:nEntries){
      nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
      snpNames[i1] <- readChar(tempCon, nChars)
    }

  }

  seek(tempCon, 15278138)
    
  normIDByteOffset <- seek(tempCon)
  normID <- readBin(tempCon, "integer", n=nEntries, size=1, signed=FALSE) + 1
  
  newBlockByteOffset <- seek(tempCon)
  newBlock <- readBin(tempCon, "integer", n=10000, size=1, signed=FALSE)
  
  close(tempCon)

  byteOffsets <- list(entriesByteOffset=entriesByteOffset,
                      #snpIndexByteOffset=snpIndexByteOffset, 
                      #snpNamesByteOffset=snpNamesByteOffset,
                      normIDByteOffset=normIDByteOffset,
                      newBlockByteOffset=newBlockByteOffset)
  
  allStuff <- list(prefixCheck=prefixCheck,
                   null.1=null.1,
                   versionNumber=versionNumber,
                   chipType=chipType,
                   null.2=null.2,
                   csvLines=csvLines,
                   nEntries=nEntries,
                   #snpIndex=snpIndex,
                   #snpNames=snpNames,
                   normID=normID,
                   newBlock=newBlock,
                   byteOffsets=byteOffsets)
  allStuff
  
}


RGtoXY = function(RG, chipType, remove=TRUE, verbose=TRUE) {
  chipList = c("human1mv1c",             # 1M
               "human370v1c",            # 370CNV
               "human650v3a",            # 650Y
               "human610quadv1b",        # 610 quad
               "human660quadv1a",        # 660 quad
               "human370quadv3c",        # 370CNV quad
               "human550v3b",            # 550K
               "human1mduov3b")          # 1M Duo   
  if(missing(chipType))
     chipType = match.arg(annotation(RG), chipList)
  else
     chipType = match.arg(chipType, chipList)
  
  pkgname <- getCrlmmAnnotationName(chipType)
  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
     suggCall <- paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
     msg <- paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
     message(strwrap(msg))
     stop("Package ", pkgname, " could not be found.")
     rm(suggCall, msg)
  }
  if(verbose) message("Loading chip annotation information.")
    loader("annotation.rda", .crlmmPkgEnv, pkgname)
#    data(annotation, package=pkgname, envir=.crlmmPkgEnv)
  annot <- getVarInEnv("annot")
  
  if(remove)
    annot = annot[annot[,"ToCall"]==1,]
  nsnps = nrow(annot)
  narrays = ncol(RG)
  aidcol = match("AddressA_ID", colnames(annot))
  if(is.na(aidcol))
    aidcol = match("Address", colnames(annot))
  bidcol = match("AddressB_ID", colnames(annot))
  if(is.na(bidcol)) 
    bidcol = match("Address2", colnames(annot))
  aids = annot[, aidcol]
  bids = annot[, bidcol]
  snpids = annot[,"Name"]
  snpbase = annot[,"SNP"] 
  rrgg = !is.na(bids) & bids!=0  
  aord = match(aids, featureNames(RG)) # NAs are possible here
  bord = match(bids[!rrgg], featureNames(RG)) # and here
#  argrg = aids[rrgg]
#  brgrg = bids[rrgg]
  X = Y = Xnb = Ynb = Xse = Yse = zero = matrix(0, nsnps, narrays)
  rownames(X) = rownames(Y) = rownames(Xnb) = rownames(Ynb) = rownames(Xse) = rownames(Yse) = snpids
  colnames(X) = colnames(Y) = colnames(Xnb) = colnames(Ynb) = colnames(Xse) = colnames(Yse) = sampleNames(RG)$G
  X[!is.na(aord),] = exprs(channel(RG, "R"))[aord[!is.na(aord)],] # mostly red
  Y[!is.na(aord),] = exprs(channel(RG, "G"))[aord[!is.na(aord)],] # mostly green
  Xnb[!is.na(aord),] = exprs(channel(RG, "Rnb"))[aord[!is.na(aord)],]
  Ynb[!is.na(aord),] = exprs(channel(RG, "Gnb"))[aord[!is.na(aord)],]
  Xse[!is.na(aord),] = exprs(channel(RG, "Rse"))[aord[!is.na(aord)],]
  Yse[!is.na(aord),] = exprs(channel(RG, "Gse"))[aord[!is.na(aord)],]
  
  # Important addition needed to get the correct intensities for the transitional Infinium I SNPs
  # these SNPs are intergorated by 2 separate bead types (the second bead type is given in the 'AddressB_ID'
  # (or 'Address2') column, and will have bid != 0 | !is.na(bid)
  # for these it could be that the X signal is the R/G from the A bead and the Y signal the R/G from the B snp
  # at present these are all zeroed out, so will presumably be called as hets
  # this is not a problem for the training data, since we get the X and Y intensities directly from BeadStudio output
  
  X[!is.na(bord),] = 0
  Y[!is.na(bord),] = 0
  Xnb[!is.na(bord),] = 0
  Ynb[!is.na(bord),] = 0
  Xse[!is.na(bord),] = 0
  Yse[!is.na(bord),] = 0

  zero[X==0 | Y==0] = 1
  
  rm(annot)
  gc()
  
  XY = new("NChannelSet",
           X=X, Y=Y,
           Xnb=Xnb, Ynb=Ynb,
           Xse=Xse, Yse=Yse,
           zero=zero,
           annotation=chipType,
           phenoData=RG@phenoData)
}

stripNormalize = function(XY, useTarget=TRUE, verbose=TRUE) {
  pkgname <- getCrlmmAnnotationName(annotation(XY))
  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
    suggCall <- paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
    msg <- paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
    message(strwrap(msg))
    stop("Package ", pkgname, " could not be found.")
    rm(suggCall, msg)
  }
  if(verbose) message("Loading strip and reference normalization information.")
  loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)
#  data(preprocStuff, package=pkgname, envir=.crlmmPkgEnv)

  stripnum <- getVarInEnv("stripnum")
  
  if(useTarget)
    targetdist = getVarInEnv("reference")
  
  Xqws = Yqws = matrix(0, nrow(XY), ncol(XY))
  colnames(Xqws) = colnames(Yqws) = sampleNames(XY)$X
  rownames(Xqws) = rownames(Yqws) = featureNames(XY)

  if(verbose){
    message("Quantile normalizing ", ncol(XY), " arrays by ", max(stripnum), " strips.")
    if (getRversion() > '2.7.0') pb <- txtProgressBar(min=0, max=max(stripnum), style=3)
  } 
  for(s in 1:max(stripnum)) {
    if(verbose) {
      if (getRversion() > '2.7.0') setTxtProgressBar(pb, s)
      else cat(".")
    }
    sel = stripnum==s
    subX = exprs(channel(XY, "X"))[sel,]
    subY = exprs(channel(XY, "Y"))[sel,]
    if(useTarget)
      tmp = normalize.quantiles.use.target(as.matrix(cbind(subX, subY)), targetdist[[s]])
    else
      tmp = normalize.quantiles(as.matrix(cbind(subX, subY)))
    Xqws[sel,] = tmp[,1:(ncol(tmp)/2)]
    Yqws[sel,] = tmp[,(ncol(tmp)/2+1):ncol(tmp)]
    rm(subX, subY, tmp, sel)
    gc()
  }
  if(verbose)
    cat("\n")
  XYNorm = new("NChannelSet",
                X=Xqws+16,
                Y=Yqws+16,
                Xnb=exprs(channel(XY, "Xnb")),
                Ynb=exprs(channel(XY, "Ynb")),
                Xse=exprs(channel(XY, "Xse")),
                Yse=exprs(channel(XY, "Yse")),
                zero=exprs(channel(XY, "zero")),
                annotation=annotation(XY),
                phenoData=XY@phenoData)
}


preprocessInfinium2 = function(XY, mixtureSampleSize=10^5, fitMixture=TRUE, eps=0.1, verbose=TRUE, seed=1,
                               cdfName, sns, stripNorm=TRUE, useTarget=TRUE) {
  if(stripNorm)
    XY = stripNormalize(XY, useTarget=useTarget, verbose=verbose)
  
## MR: the code below is mostly straight from snprma.R
  if (missing(sns)) sns <- sampleNames(XY)
  if(missing(cdfName))
    cdfName <- annotation(XY)
##  stuffDir <- changeToCrlmmAnnotationName(cdfName)
  pkgname <- getCrlmmAnnotationName(cdfName)
  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
    suggCall <- paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
    msg <- paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
    message(strwrap(msg))
    stop("Package ", pkgname, " could not be found.")
    rm(suggCall, msg)
  }
  if(verbose) message("Loading snp annotation and mixture model parameters.")
  loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
  loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)  
#  data(genotypeStuff, mixtureStuff, package=pkgname, envir=.crlmmPkgEnv)
  autosomeIndex <- getVarInEnv("autosomeIndex")
  #  pnsa <- getVarInEnv("pnsa")
  #  pnsb <- getVarInEnv("pnsb")
  #  fid <- getVarInEnv("fid")
  #  reference <- getVarInEnv("reference")
  #  aIndex <- getVarInEnv("aIndex")
  #  bIndex <- getVarInEnv("bIndex")
  SMEDIAN <- getVarInEnv("SMEDIAN")
  theKnots <- getVarInEnv("theKnots")
  gns <- getVarInEnv("gns")

  nprobes <- nrow(XY)
  narrays <- ncol(XY)
  
  ##We will read each cel file, summarize, and run EM one by one
  ##We will save parameters of EM to use later
  mixtureParams <- matrix(0, 4, narrays)
  SNR <- vector("numeric", narrays)
  SKW <- vector("numeric", narrays)

  ## This is the sample for the fitting of splines
  ## BC: I like better the idea of the user passing the seed,
  ##     because this might intefere with other analyses
  ##     (like what happened to GCRMA)
  set.seed(seed)

  idx <- sort(sample(autosomeIndex, mixtureSampleSize))
  
  ##S will hold (A+B)/2 and M will hold A-B
  ##NOTE: We actually dont need to save S. Only for pics etc...
  ##f is the correction. we save to avoid recomputing
  A <- matrix(as.integer(exprs(channel(XY, "X"))), nprobes, narrays) # matrix(as.integer(0), length(pnsa), length(filenames))
  B <- matrix(as.integer(exprs(channel(XY, "Y"))), nprobes, narrays) # matrix(as.integer(0), length(pnsb), length(filenames))

  # new lines below - useful to keep track of zeroed out SNPs
  zero <- matrix(as.integer(exprs(channel(XY, "zero"))), nprobes, narrays)
    
  colnames(A) <- colnames(B) <- colnames(zero) <- sampleNames(XY)$X
  rownames(A) <- rownames(B) <- rownames(zero) <- featureNames(XY)
  
  if(verbose){
     message("Calibrating ", narrays, " arrays.")
     if (getRversion() > '2.7.0') pb <- txtProgressBar(min=0, max=narrays, style=3)
  }

  idx2 <- sample(nprobes, 10^5)
  for(i in 1:narrays){
    SKW[i] = mean((A[idx2,i]-mean(A[idx2,i]))^3)/(sd(A[idx2,i])^3)
    if(fitMixture){
      S <- (log2(A[idx, i])+log2(B[idx, i]))/2 - SMEDIAN
      M <- log2(A[idx, i])-log2(B[idx, i])

      ##we need to test the choice of eps.. it is not the max diff between funcs
      tmp <- fitAffySnpMixture56(S, M, theKnots, eps=eps)

      mixtureParams[, i] <- tmp[["coef"]]
      SNR[i] <- tmp[["medF1"]]^2/(tmp[["sigma1"]]^2+tmp[["sigma2"]]^2)
    }
    if(verbose) {
      if (getRversion() > '2.7.0') setTxtProgressBar(pb, i)
      else cat(".")
    }
  }
  if (verbose) {
    if (getRversion() > '2.7.0') close(pb)
    else cat("\n")
  }
  if (!fitMixture) SNR <- mixtureParams <- NA
  ## gns comes from preprocStuff.rda
  list(A=A, B=B, zero=zero, sns=sns, gns=gns, SNR=SNR, SKW=SKW, mixtureParams=mixtureParams, cdfName=cdfName)
}


## MR: Could add arguments to allow this function to read in idat data as well, although this add a further 7 arguments, which might over-complicate things
crlmmIllumina <- function(RG, XY, stripNorm=TRUE, useTarget=TRUE,
                  row.names=TRUE, col.names=TRUE,
                  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                  seed=1, # save.it=FALSE, load.it=FALSE, intensityFile,
                  mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
                  cdfName, sns, recallMin=10, recallRegMin=1000,
                  returnParams=FALSE, badSNP=.7) {
  if(!missing(RG)) {
    if(missing(XY))
      XY = RGtoXY(RG, chipType=cdfName)
    else
      stop("Both RG and XY specified - please use one or the other")
  }
#  if ((load.it | save.it) & missing(intensityFile))
#    stop("'intensityFile' is missing, and you chose either load.it or save.it")
  if (missing(sns)) sns <- sampleNames(XY) #basename(filenames)
#  if (!missing(intensityFile))
#    if (load.it & !file.exists(intensityFile)){
#      load.it <- FALSE
#      message("File ", intensityFile, " does not exist.")
#      message("Not loading it, but running SNPRMA from scratch.")
# }
#  if (!load.it){
#    res <- snprma(filenames, fitMixture=TRUE,
#                    mixtureSampleSize=mixtureSampleSize, verbose=verbose,
#                    eps=eps, cdfName=cdfName, sns=sns)
  res = preprocessInfinium2(XY, mixtureSampleSize=mixtureSampleSize, fitMixture=TRUE, verbose=verbose,
                        seed=seed, eps=eps, cdfName=cdfName, sns=sns, stripNorm=stripNorm, useTarget=useTarget)
#    if(save.it){
#      t0 <- proc.time()
#      save(res, file=intensityFile)
#      t0 <- proc.time()-t0
#      if (verbose) message("Used ", t0[3], " seconds to save ", intensityFile, ".")
#    }
#      }else{
#            if (verbose) message("Loading ", intensityFile, ".")
#                obj <- load(intensityFile)
#                if (verbose) message("Done.")
#                if (obj != "res")
#                        stop("Object in ", intensityFile, " seems to be invalid.")
#          }
  if(row.names) row.names=res$gns else row.names=NULL
  if(col.names) col.names=res$sns else col.names=NULL

  res2 <- crlmmGT(res[["A"]], res[["B"]], res[["SNR"]],
                  res[["mixtureParams"]], res[["cdfName"]],
                  gender=gender, row.names=row.names,
                  col.names=col.names, recallMin=recallMin,
                  recallRegMin=1000, SNRMin=SNRMin,
                  returnParams=returnParams, badSNP=badSNP,
                  verbose=verbose)
  
  res2[["A"]] <- res[["A"]]  # added for copy number analysis
  res2[["B"]] <- res[["B"]]  # added for copy number analysis
  res2[["SNR"]] <- res[["SNR"]]
  res2[["SKW"]] <- res[["SKW"]]
  res2[["zero"]] <- res[["zero"]]
  rm(res)
  gc()
  # MR: FIXME - for consistency with crlmm, need to save results in a 'SnpSet' object
  return(res2)
}
