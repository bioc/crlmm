# function below works OK provided all .idat files are in the current working directory
# - could add an option to allow files in Illumina directory structure to be handled
# or to use the optional 'Path' column in sampleSheet
# - there is a lot of header information that is currently discarded - could try and store this somewhere in the resulting NChannelSet

readIdatFiles = function(sampleSheet=NULL,
			  arrayNames=NULL,
			  ids=NULL,
			  path=".",
			  arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			  highDensity=FALSE,
			  sep="_",
			  fileExt=list(green="Grn.idat", red="Red.idat"),
			  saveDate=FALSE) {
       if(!is.null(arrayNames)) {
               pd = new("AnnotatedDataFrame", data = data.frame(Sample_ID=arrayNames))
       }
       if(!is.null(sampleSheet)) { # get array info from Illumina's sample sheet
	       if(is.null(arrayNames)){
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
	       }
	       pd = new("AnnotatedDataFrame", data = sampleSheet)
	       sampleNames(pd) = basename(arrayNames)
       }
       if(is.null(arrayNames)) {
               arrayNames = gsub(paste(sep, fileExt$green, sep=""), "", dir(pattern=fileExt$green, path=path))
               if(!is.null(sampleSheet)) {
                      sampleSheet=NULL
                      cat("Could not find required info in \'sampleSheet\' - ignoring.  Check \'sampleSheet\' and/or \'arrayInfoColNames\'\n")
               }
               pd = new("AnnotatedDataFrame", data = data.frame(Sample_ID=arrayNames))
       }

       narrays = length(arrayNames)
       grnfiles = paste(arrayNames, fileExt$green, sep=sep)
       redfiles = paste(arrayNames, fileExt$red, sep=sep)
       if(length(grnfiles)==0 || length(redfiles)==0)
	       stop("Cannot find .idat files")
       if(length(grnfiles)!=length(redfiles))
	       stop("Cannot find matching .idat files")
       if(path[1] != "." & path[1] != ""){
	       grnidats = file.path(path, grnfiles)
	       redidats = file.path(path, redfiles)
       }  else {
	       message("path arg not set.  Assuming files are in local directory, or that complete path is provided")
	       grnidats = grnfiles
	       redidats = redfiles
       }
       if(!all(file.exists(grnidats))) stop("Missing some of the *Grn.idat files")
       if(!all(file.exists(redidats))) stop("Missing some of the *Red.idat files")
       headerInfo = list(nProbes = rep(NA, narrays),
                         Barcode = rep(NA, narrays),
                         ChipType = rep(NA, narrays),
                         Manifest = rep(NA, narrays), # not sure about this one - sometimes blank
                         Position = rep(NA, narrays)) # this may also vary a bit
       dates = list(decode=rep(NA, narrays),
                    scan=rep(NA, narrays))
       ## read in the data
       for(i in seq_along(arrayNames)) {
	       cat("reading", arrayNames[i], "\t")
	       idsG = idsR = G = R = NULL
	       cat(paste(sep, fileExt$green, sep=""), "\t")
	       G = readIDAT(grnidats[i])
	       idsG = rownames(G$Quants)
	       headerInfo$nProbes[i] = G$nSNPsRead
	       headerInfo$Barcode[i] = G$Barcode
	       headerInfo$ChipType[i] = G$ChipType
	       headerInfo$Manifest[i] = G$Unknown$MostlyNull
	       headerInfo$Position[i] = G$Unknowns$MostlyA
               if(headerInfo$nProbes[i]>(headerInfo$nProbes[1]+10000) || headerInfo$nProbes[i]<(headerInfo$nProbes[1]-10000)) {
		       warning("Chips are not of the same type.  Skipping ", basename(grnidats[i]), " and ", basename(redidats[i]))
		       next()
	       }
	       dates$decode[i] = G$RunInfo[1, 1]
	       dates$scan[i] = G$RunInfo[2, 1]
	       if(i==1) {
		       if(is.null(ids) && !is.null(G)){
			       ids = idsG
		       } else stop("Could not find probe IDs")
		       nprobes = length(ids)
		       narrays = length(arrayNames)
		       RG = new("NChannelSet",
		                 R=initializeBigMatrix(name="R", nr=nprobes, nc=narrays, vmode="integer"),
		                 G=initializeBigMatrix(name="G", nr=nprobes, nc=narrays, vmode="integer"),
		                 zero=initializeBigMatrix(name="zero", nr=nprobes, nc=narrays, vmode="integer"),
				         annotation=headerInfo$Manifest[1],
				         phenoData=pd, storage.mode="environment")
			   featureNames(RG) = ids
		       if(!is.null(sampleSheet) && !is.null(sampleSheet$Sample_ID)){
		            sampleNames(RG) = sampleSheet$Sample_ID
		       } else  sampleNames(RG) = arrayNames
		       gc()
	       }
	       if(length(ids)==length(idsG)) {
		       if(sum(ids==idsG)==nprobes) {
			       RG@assayData$G[,i] = G$Quants[, "Mean"]
			       zeroG = G$Quants[, "NBeads"]==0
		       }
	       } else {
		       indG = match(ids, idsG)
		       RG@assayData$G[,i] = G$Quants[indG, "Mean"]
		       zeroG = G$Quants[indG, "NBeads"]==0
	       }
	       rm(G)
	       gc()

	       cat(paste(sep, fileExt$red, sep=""), "\n")
	       R = readIDAT(redidats[i])
	       idsR = rownames(R$Quants)

	       if(length(ids)==length(idsG)) {
		       if(sum(ids==idsR)==nprobes) {
			       RG@assayData$R[,i] = R$Quants[ ,"Mean"]
		           zeroR = R$Quants[ ,"NBeads"]==0
		       }
	       } else {
		       indR = match(ids, idsR)
		       RG@assayData$R[,i] = R$Quants[indR, "Mean"]
		       zeroR = R$Quants[indR, "NBeads"]==0
	       }
	       RG@assayData$zero[,i] = zeroG | zeroR
	       rm(R, zeroG, zeroR)
	       gc()
       }
       if(saveDate) {
	       protocolData(RG)[["ScanDate"]] = dates$scan
       }
       storageMode(RG) = "lockedEnvironment"
       is.lds = ifelse(isPackageLoaded("ff"), TRUE, FALSE)       
       if(is.lds) {
         close(RG@assayData$R)
         close(RG@assayData$G)
         close(RG@assayData$zero)
       }
       RG
}


## the readIDAT() and readBPM() functions below were provided by Keith Baggerly, 27/8/2008
readIDAT <- function(idatFile){
  fileSize <- file.info(idatFile)$size

  tempCon <- file(idatFile,"rb")
  prefixCheck <- readChar(tempCon,4)
  if(prefixCheck != "IDAT"){

  }

  versionNumber <- readBin(tempCon, "integer", n=1, size=8,
                           endian="little", signed=FALSE)

  if(versionNumber<3)
	  stop("Older style IDAT files not supported:  consider updating your scanner settings")

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
  for(i1 in 1:2) { #nRunInfoBlocks){  ## MR edit
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


RGtoXY = function(RG, chipType, verbose=TRUE) {
  chipList = c("human1mv1c",             # 1M
               "human370v1c",            # 370CNV
               "human650v3a",            # 650Y
               "human610quadv1b",        # 610 quad
               "human660quadv1a",        # 660 quad
               "human370quadv3c",        # 370CNV quad
               "human550v3b",            # 550K
               "human1mduov3b",          # 1M Duo
               "humanomni1quadv1b",      # Omni1 quad
	       "humanomniexpress12v1b",  # Omni express 12
               "humanimmuno12v1b")       # Immuno chip 12
  if(missing(chipType)){
	  chipType = match.arg(annotation(RG), chipList)
  } else chipType = match.arg(chipType, chipList)

  pkgname = getCrlmmAnnotationName(chipType)
  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
     suggCall = paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
     msg = paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
     message(strwrap(msg))
     stop("Package ", pkgname, " could not be found.")
     rm(suggCall, msg)
  }
  if(verbose) message("Loading chip annotation information.")
    loader("address.rda", .crlmmPkgEnv, pkgname)
  aids = getVarInEnv("addressA") # comes from AddressA_ID or Address column in manifest
  bids = getVarInEnv("addressB") # comes from AddressB_ID or Address2 column in manifest
  ids = names(aids)
  snpbase = getVarInEnv("base")

  nsnps = length(aids)
  narrays = ncol(RG)

#  aidcol = match("AddressA_ID", colnames(annot))
#  if(is.na(aidcol))
#    aidcol = match("Address", colnames(annot))
#  bidcol = match("AddressB_ID", colnames(annot))
#  if(is.na(bidcol))
#    bidcol = match("Address2", colnames(annot))
#  aids = annot[, aidcol]
#  bids = annot[, bidcol]
#  snpids = annot[,"Name"]
#  snpbase = annot[,"SNP"]
  infI = !is.na(bids) & bids!=0
  aord = match(aids, featureNames(RG)) # NAs are possible here
  bord = match(bids, featureNames(RG)) # and here
#  argrg = aids[rrgg]
#  brgrg = bids[rrgg]

  XY = new("NChannelSet",
	     X=initializeBigMatrix(name="X", nr=nsnps, nc=narrays, vmode="integer"),
	     Y=initializeBigMatrix(name="Y", nr=nsnps, nc=narrays, vmode="integer"),
	     zero=initializeBigMatrix(name="zero", nr=nsnps, nc=narrays, vmode="integer"),
	     annotation=chipType, phenoData=RG@phenoData,
	     protocolData=RG@protocolData, storage.mode="environment")
  featureNames(XY) = ids 
  sampleNames(XY) = sampleNames(RG)
  gc()
  # Need to initialize - matrices filled with NAs to begin with
  XY@assayData$X[1:nsnps,] = 0
  XY@assayData$Y[1:nsnps,] = 0
  XY@assayData$zero[1:nsnps,] = 0

  is.lds = ifelse(isPackageLoaded("ff"), TRUE, FALSE)
  
  if(is.lds) {
    open(RG@assayData$G)
    open(RG@assayData$R)
    open(RG@assayData$zero)
  }

  # First sort out Infinium II SNPs, X -> R (allele A)  and Y -> G (allele B) from the same probe
  XY@assayData$X[!is.na(aord),] = exprs(channel(RG, "R"))[aord[!is.na(aord)],] # mostly red
  XY@assayData$Y[!is.na(aord),] = exprs(channel(RG, "G"))[aord[!is.na(aord)],] # mostly green
  XY@assayData$zero[!is.na(aord),] = exprs(channel(RG, "zero"))[aord[!is.na(aord)],] # mostly green
  gc()

  ## Warning - not 100% sure that the code below is correct - could be more complicated than this

  # Next Infinium I where X -> R from allele A probe and Y -> R from allele B probe
#  infIRR = infI & (snpbase=="[A/T]" | snpbase=="[T/A]" | snpbase=="[a/t]" | snpbase=="[t/a]")

#  X[infIRR,] = exprs(channel(RG, "R"))[aord[infIRR],] # mostly red
#  Y[infIRR,] = exprs(channel(RG, "R"))[bord[infIRR],] # mostly green

  # Finally Infinium I where X -> G from allele A probe and Y -> G from allele B probe
#  infIGG = infI & (snpbase=="[C/G]" | snpbase=="[G/C]" | snpbase=="[g/c]" | snpbase=="[c/g]")

#  X[infIGG,] = exprs(channel(RG, "G"))[aord[infIGG],] # mostly red
#  Y[infIGG,] = exprs(channel(RG, "G"))[bord[infIGG],] # mostly green

  #  For now zero out Infinium I probes
  XY@assayData$X[infI,] = 0
  XY@assayData$Y[infI,] = 0
  XY@assayData$zero[infI,] = 0
  gc()

  is.lds = ifelse(isPackageLoaded("ff"), TRUE, FALSE)
  if(is.lds) {
    close(RG@assayData$G)
    close(RG@assayData$R)
    close(RG@assayData$zero)    
    close(XY@assayData$X)
    close(XY@assayData$Y)
    close(XY@assayData$zero)
  }
  XY
}


stripNormalize = function(XY, useTarget=TRUE, verbose=TRUE) {
  pkgname = getCrlmmAnnotationName(annotation(XY))
  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
    suggCall = paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
    msg = paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
    message(strwrap(msg))
    stop("Package ", pkgname, " could not be found.")
    rm(suggCall, msg)
  }
  if(verbose) message("Loading strip and reference normalization information.")
  loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)

  stripnum = getVarInEnv("stripnum")

  if(useTarget)
    targetdist = getVarInEnv("reference")

  if(verbose){
    message("Quantile normalizing ", ncol(XY), " arrays by ", max(stripnum), " strips.")
    if (getRversion() > '2.7.0') pb = txtProgressBar(min=0, max=max(stripnum), style=3)
  }
  is.lds = ifelse(isPackageLoaded("ff"), TRUE, FALSE)
  
  if(is.lds) {
    open(XY@assayData$X)
    open(XY@assayData$Y)
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
    XY@assayData$X[sel,] = matrix(as.integer(tmp[,1:(ncol(tmp)/2)]+16), nrow(tmp), ncol(tmp)/2)
    XY@assayData$Y[sel,] = matrix(as.integer(tmp[,(ncol(tmp)/2+1):ncol(tmp)]+16), nrow(tmp), ncol(tmp)/2)
    rm(subX, subY, tmp, sel)
    gc()
  }
  if(is.lds) {
    close(XY@assayData$X)
    close(XY@assayData$Y)
  }

  if(verbose)
    cat("\n")
  XY
}


preprocessInfinium2 = function(XY, mixtureSampleSize=10^5,
				fitMixture=TRUE,
				eps=0.1,
				verbose=TRUE,
				seed=1,
				cdfName,
				sns,
				stripNorm=TRUE,
				useTarget=TRUE) { 
#				save.it=FALSE,
#				snpFile,
#				cnFile) {
  if(stripNorm)
    XY = stripNormalize(XY, useTarget=useTarget, verbose=verbose)

## MR: the code below is mostly straight from snprma.R
  if (missing(sns)) sns = sampleNames(XY) #$X
  if(missing(cdfName))
    cdfName = annotation(XY)
  pkgname = getCrlmmAnnotationName(cdfName)
  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
    suggCall = paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
    msg = paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
    message(strwrap(msg))
    stop("Package ", pkgname, " could not be found.")
    rm(suggCall, msg)
  }
  if(verbose) message("Loading snp annotation and mixture model parameters.")
  loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
  loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)
  loader("snpProbesFid.rda", .crlmmPkgEnv, pkgname)
  loader("npProbesFid.rda", .crlmmPkgEnv, pkgname)
  autosomeIndex = getVarInEnv("autosomeIndex")
  SMEDIAN = getVarInEnv("SMEDIAN")
  theKnots = getVarInEnv("theKnots")
  narrays = ncol(XY)

  is.lds = ifelse(isPackageLoaded("ff"), TRUE, FALSE)    
  
#  if(save.it & !missing(cnFile)) {
    # separate out copy number probes
    npIndex = getVarInEnv("npProbesFid")
    nprobes = length(npIndex)
    if(length(nprobes)>0) {
      if(is.lds) {    
        open(XY@assayData$X)
        open(XY@assayData$Y)
        open(XY@assayData$zero)
      }
      A = matrix(as.integer(exprs(channel(XY, "X"))[npIndex,]), nprobes, narrays)
      B = matrix(as.integer(exprs(channel(XY, "Y"))[npIndex,]), nprobes, narrays)

      # new lines below - useful to keep track of zeroed out probes
      zero = matrix(as.integer(exprs(channel(XY, "zero"))[npIndex,]), nprobes, narrays)

      colnames(A) = colnames(B) = colnames(zero) = sns
      rownames(A) = rownames(B) = rownames(zero) = names(npIndex)

      cnAB = list(A=A, B=B, zero=zero, sns=sns, gns=names(npIndex), cdfName=cdfName)

#      t0 <- proc.time()
#      save(cnAB, file=cnFile)
#      t0 <- proc.time()-t0
#      if(verbose) message("Used ", round(t0[3],1), " seconds to save ", cnFile, ".")
       rm(A, B, zero)
      if(is.lds) {    
        close(XY@assayData$X)
        close(XY@assayData$Y)
        close(XY@assayData$zero)
      }
    }
#  }

  # next process snp probes
  snpIndex = getVarInEnv("snpProbesFid")
  nprobes = length(snpIndex)

  ##We will read each cel file, summarize, and run EM one by one
  ##We will save parameters of EM to use later
  mixtureParams = initializeBigMatrix("crlmmMixt-", 4, narrays, "double")
  SNR = initializeBigVector("crlmmSNR-", narrays, "double")
  SKW = initializeBigVector("crlmmSKW-", narrays, "double")

  ## This is the sample for the fitting of splines
  ## BC: I like better the idea of the user passing the seed,
  ##     because this might intefere with other analyses
  ##     (like what happened to GCRMA)
  set.seed(seed)
  idx = sort(sample(autosomeIndex, mixtureSampleSize))
  idx2 = sample(nprobes, 10^5)

  ##S will hold (A+B)/2 and M will hold A-B
  ##NOTE: We actually dont need to save S. Only for pics etc...
  ##f is the correction. we save to avoid recomputing

  A = initializeBigMatrix("crlmmA-", nprobes, narrays, "integer")
  B = initializeBigMatrix("crlmmB-", nprobes, narrays, "integer")
  zero = initializeBigMatrix("crlmmZero-", nprobes, narrays, "integer")

  if(verbose){
     message("Calibrating ", narrays, " arrays.")
     if (getRversion() > '2.7.0') pb = txtProgressBar(min=0, max=narrays, style=3)
  }

  if(is.lds) {
    open(XY@assayData$X)
    open(XY@assayData$Y)
    open(XY@assayData$zero)
  }
  
  for(i in 1:narrays){
     A[,i] = as.integer(exprs(channel(XY, "X"))[snpIndex,i])
     B[,i] = as.integer(exprs(channel(XY, "Y"))[snpIndex,i])
     zero[,i] = as.integer(exprs(channel(XY, "zero"))[snpIndex,i])
     SKW[i] = mean((A[idx2,i]-mean(A[idx2,i]))^3)/(sd(A[idx2,i])^3)
     if(fitMixture){
        S = (log2(A[idx,i])+log2(B[idx,i]))/2 - SMEDIAN
        M = log2(A[idx,i])-log2(B[idx,i])

        ##we need to test the choice of eps.. it is not the max diff between funcs
        tmp = fitAffySnpMixture56(S, M, theKnots, eps=eps)

        mixtureParams[, i] = tmp[["coef"]]
        SNR[i] = tmp[["medF1"]]^2/(tmp[["sigma1"]]^2+tmp[["sigma2"]]^2)
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
  if (!fitMixture) SNR = mixtureParams = NA
  ## gns comes from preprocStuff.rda

#  if(save.it & !missing(snpFile)) {
#    t0 <- proc.time()
#    save(res, file=snpFile)
#    t0 <- proc.time()-t0
#    if(verbose) message("Used ", round(t0[3],1), " seconds to save ", snpFile, ".")
#  }
  if(is.lds) {
    close(XY@assayData$X)
    close(XY@assayData$Y)
    close(XY@assayData$zero)    
    close(A)
    close(B)
    close(zero)
    close(SKW)
    close(mixtureParams)
    close(SNR)
  }

  res = list(A=A, B=B,
             zero=zero, sns=sns, gns=names(snpIndex), SNR=SNR, SKW=SKW,
             mixtureParams=mixtureParams, cdfName=cdfName, cnAB=cnAB)
  
#  open(res[["A"]])
#  open(res[["B"]])
#  open(res[["zero"]])
#  open(res[["SNR"]])
#  open(res[["mixtureParams"]])

#  if(save.it & !missing(snpFile)) {
#    t0 <- proc.time()
#    save(res, file=snpFile)
#    t0 <- proc.time()-t0
#    if(verbose) message("Used ", round(t0[3],1), " seconds to save ", snpFile, ".")
#  }

#  close(res[["A"]])
#  close(res[["B"]])
#  close(res[["zero"]])
#  close(res[["SNR"]])
#  close(res[["mixtureParams"]])

  return(res)
}


crlmmIllumina = function(RG, XY, stripNorm=TRUE, useTarget=TRUE,
                  row.names=TRUE, col.names=TRUE,
                  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                  seed=1, # save.it=FALSE, load.it=FALSE, snpFile, cnFile,
                  mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
                  cdfName, sns, recallMin=10, recallRegMin=1000,
                  returnParams=FALSE, badSNP=.7) {
#  if (save.it & (missing(snpFile) | missing(cnFile)))
#    stop("'snpFile' and/or 'cnFile' is missing and you chose to save.it")
#  if (load.it & missing(snpFile))
#    stop("'snpFile' is missing and you chose to load.it")
#  if (!missing(snpFile))
#    if (load.it & !file.exists(snpFile)){
#      load.it <- FALSE
#      message("File ", snpFile, " does not exist.")
#      stop("Cannot load SNP data.")
#  }
#  if (!load.it){
    if(!missing(RG)) {
      if(missing(XY))
        XY = RGtoXY(RG, chipType=cdfName)
      else
        stop("Both RG and XY specified - please use one or the other")
    }
    if (missing(sns)) sns = sampleNames(XY)
    res = preprocessInfinium2(XY, mixtureSampleSize=mixtureSampleSize, fitMixture=TRUE, verbose=verbose,
                        seed=seed, eps=eps, cdfName=cdfName, sns=sns, stripNorm=stripNorm, useTarget=useTarget) #,
#                        save.it=save.it, snpFile=snpFile, cnFile=cnFile)

    is.lds = ifelse(isPackageLoaded("ff"), TRUE, FALSE)

    if(is.lds) {
      open(res[["A"]])
      open(res[["B"]])
      open(res[["SNR"]])
      open(res[["mixtureParams"]])
    }
    
#    fD = featureData(XY)
#    phenD = XY@phenoData
#    protD = XY@protocolData
#    rm(XY)
#    gc()
#    if(verbose) message("Initializing container for alleleA, alleleB, call, callProbability")
#    callSet <- new("SnpSuperSet",
#                    alleleA=initializeBigMatrix(name="A", nr=nrow(res[[1]]), nc=length(sns)),
#                    alleleB=initializeBigMatrix(name="B", nr=nrow(res[[1]]), nc=length(sns)),
#                    call=initializeBigMatrix(name="call", nr=nrow(res[[1]]), nc=length(sns)),
#                    callProbability=initializeBigMatrix(name="callPr", nr=nrow(res[[1]]), nc=length(sns)),
#  	            annotation=cdfName, protocolData=protD, phenoData=phenD, featureData=fD)
#    sampleNames(callSet) <- sns
#    featureNames(callSet) <- res[["gns"]]
#    pData(callSet)$SKW <- rep(NA, length(sns))
#    pData(callSet)$SNR <- rep(NA, length(sns))
#    pData(callSet)$gender <- rep(NA, length(sns))

#  }else{
#      if(verbose) message("Loading ", snpFile, ".")
#        obj <- load(snpFile)
#        if(verbose) message("Done.")
#        if(!any(obj == "res"))
#          stop("Object in ", snpFile, " seems to be invalid.")
#  }

 #   rm(phenD, protD , fD)

#    snp.index <- res$snpIndex #match(res$gns, featureNames(callSet))
#    suppressWarnings(A(callSet) <- res[["A"]])
#    suppressWarnings(B(callSet) <- res[["B"]])
#    pData(callSet)$SKW <- res$SKW
#    pData(callSet)$SNR <- res$SNR
#    mixtureParams <- res$mixtureParams
#    rm(res); gc()
  if(row.names) row.names=res$gns else row.names=NULL
  if(col.names) col.names=res$sns else col.names=NULL

  FUN = ifelse(is.lds, "crlmmGT2", "crlmmGT")
  ## genotyping
  crlmmGTfxn = function(FUN,...){
		switch(FUN,
		       crlmmGT2=crlmmGT2(...),
		       crlmmGT=crlmmGT(...))
              }
  res2 = crlmmGTfxn(FUN,
                     A=res[["A"]],
                     B=res[["B"]],
                     SNR=res[["SNR"]],
                     mixtureParams=res[["mixtureParams"]],
                     cdfName=cdfName,
                     row.names=row.names,
                     col.names=col.names,
                     probs=probs,
                     DF=DF,
                     SNRMin=SNRMin,
                     recallMin=recallMin,
                     recallRegMin=recallRegMin,
                     gender=gender,
                     verbose=verbose,
                     returnParams=returnParams,
                     badSNP=badSNP)
    
#  res2 <- crlmmGT2(A=res[["A"]], #as.matrix(A(callSet)[snp.index,]), # j]),
#                  B=res[["B"]], # as.matrix(B(callSet)[snp.index,]),  # j]),
#                  SNR=res[["SNR"]], # callSet$SNR, # [j],
#                  mixtureParams=res[["mixtureParams"]], #,
#                  cdfName=res[["cdfName"]], # annotation(callSet),
#                  row.names=row.names, # featureNames(callSet)[snp.index],
#                  col.names=col.names, # sampleNames(callSet), #[j],
#                  probs=probs,
#                  DF=DF,
#                  SNRMin=SNRMin,
#                  recallMin=recallMin,
#                  recallRegMin=recallRegMin,
#                  gender=gender,
#                  verbose=verbose,
#                  returnParams=returnParams,
#                  badSNP=badSNP)
#    rm(res); gc()
#    suppressWarnings(snpCall(callSet)[snp.index, j] <- tmp[["calls"]])
#    suppressWarnings(snpCallProbability(callSet)[snp.index, j] <- tmp[["confs"]])
#    callSet$gender[j] <- tmp$gender
#    suppressWarnings(snpCall(callSet)[snp.index,] <- tmp[["calls"]])
#    suppressWarnings(snpCallProbability(callSet)[snp.index,] <- tmp[["confs"]])
#    callSet$gender <- tmp$gender
#    rm(tmp); gc()
#    return(callSet)

  res2[["SNR"]] = res[["SNR"]]
  res2[["SKW"]] = res[["SKW"]]
#  if(is.lds) {
#    delete(res[["A"]])
#    delete(res[["B"]])
#    delete(res[["SNR"]])
#    delete(res[["mixtureParams"]])
#  }
  rm(res); gc()
  return(list2SnpSet(res2, returnParams=returnParams))
}


## MR: Below is a more memory efficient version of crlmmIllumina() which
## reads in the .idats and genotypes in the one function and removes objects
## after they have been used
crlmmIlluminaV2 = function(sampleSheet=NULL,
			  arrayNames=NULL,
			  ids=NULL,
			  path=".",
			  arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			  highDensity=FALSE,
			  sep="_",
			  fileExt=list(green="Grn.idat", red="Red.idat"),
			  saveDate=FALSE,
			  stripNorm=TRUE,
			  useTarget=TRUE,
			  row.names=TRUE,
			  col.names=TRUE,
			  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                          seed=1,  # save.it=FALSE, snpFile, cnFile,
                          mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
                          cdfName, sns, recallMin=10, recallRegMin=1000,
                          returnParams=FALSE, badSNP=.7) {

    if(missing(cdfName)) stop("must specify cdfName")
    if(!isValidCdfName(cdfName)) stop("cdfName not valid.  see validCdfNames")
  
    is.lds = ifelse(isPackageLoaded("ff"), TRUE, FALSE)

    RG = readIdatFiles(sampleSheet=sampleSheet, arrayNames=arrayNames,
                       ids=ids, path=path, arrayInfoColNames=arrayInfoColNames,
                       highDensity=highDensity, sep=sep, fileExt=fileExt, saveDate=saveDate)


    XY = RGtoXY(RG, chipType=cdfName)
    if(is.lds) {
      open(RG@assayData$R); open(RG@assayData$G); open(RG@assayData$zero)
      delete(RG@assayData$R); delete(RG@assayData$G); delete(RG@assayData$zero)
    }
    rm(RG); gc()
  
    if (missing(sns)) { sns = sampleNames(XY)
    } 
    res = preprocessInfinium2(XY, mixtureSampleSize=mixtureSampleSize, fitMixture=TRUE, verbose=verbose,
                               seed=seed, eps=eps, cdfName=cdfName, sns=sns, stripNorm=stripNorm, useTarget=useTarget) #,
#                               save.it=save.it, snpFile=snpFile, cnFile=cnFile)

    if(is.lds) {
      open(XY@assayData$X); open(XY@assayData$Y); open(XY@assayData$zero)
      delete(XY@assayData$X); delete(XY@assayData$Y); delete(XY@assayData$zero)
    }
    rm(XY); gc()

    if(row.names) row.names=res$gns else row.names=NULL
    if(col.names) col.names=res$sns else col.names=NULL
  
    FUN = ifelse(is.lds, "crlmmGT2", "crlmmGT")
    ## genotyping
    crlmmGTfxn = function(FUN,...){
		switch(FUN,
		       crlmmGT2=crlmmGT2(...),
		       crlmmGT=crlmmGT(...))
              }
    res2 = crlmmGTfxn(FUN,
                     A=res[["A"]],
                     B=res[["B"]],
                     SNR=res[["SNR"]],
                     mixtureParams=res[["mixtureParams"]],
                     cdfName=cdfName,
                     row.names=row.names,
                     col.names=col.names,
                     probs=probs,
                     DF=DF,
                     SNRMin=SNRMin,
                     recallMin=recallMin,
                     recallRegMin=recallRegMin,
                     gender=gender,
                     verbose=verbose,
                     returnParams=returnParams,
                     badSNP=badSNP)
  
    if(is.lds) {
      open(res[["SNR"]]); open(res[["SKW"]])
    }
    res2[["SNR"]] = res[["SNR"]]
    res2[["SKW"]] = res[["SKW"]]
 #  if(is.lds) {
 #    delete(res[["A"]]); delete(res[["B"]])
 #    delete(res[["SKW"]]); delete(res[["SNR"]]); delete(res[["mixtureParams"]])
 #  }
    rm(res); gc()
    return(list2SnpSet(res2, returnParams=returnParams))
}

# Functions analogous to Rob's Affy functions to set up container
getProtocolData.Illumina = function(filenames, sep="_", fileExt="Grn.idat") {
       narrays = length(filenames)

       headerInfo = list(nProbes = rep(NA, narrays),
                         Barcode = rep(NA, narrays),
                         ChipType = rep(NA, narrays),
                         Manifest = rep(NA, narrays),
                         Position = rep(NA, narrays))

       scanDates = data.frame(ScanDate=rep(NA, narrays), DecodeDate=rep(NA, narrays))
       rownames(scanDates) = gsub(paste(sep, fileExt, sep=""), "", filenames)
       ## read in the data
       for(i in seq_along(filenames)) {
	       cat("reading", filenames[i], "\n")
	       idsG = G = NULL
	       G = readIDAT(filenames[i])
	       idsG = rownames(G$Quants)
	       headerInfo$nProbes[i] = G$nSNPsRead
	       headerInfo$Barcode[i] = G$Barcode
	       headerInfo$ChipType[i] = G$ChipType
	       headerInfo$Manifest[i] = G$Unknown$MostlyNull
	       headerInfo$Position[i] = G$Unknowns$MostlyA
               if(headerInfo$nProbes[i]>(headerInfo$nProbes[1]+10000) || headerInfo$nProbes[i]<(headerInfo$nProbes[1]-10000)) {
		       warning("Chips are not of the same type.  Skipping ", basename(filenames[i]))
		       next()
	       }
	       scanDates$ScanDate[i] = G$RunInfo[1, 1]
	       scanDates$DecodeDate[i] = G$RunInfo[2, 1]
	       rm(G)
	       gc()
       }
       protocoldata = new("AnnotatedDataFrame",
			    data=scanDates,
			    varMetadata=data.frame(labelDescription=colnames(scanDates),
			                           row.names=colnames(scanDates)))
	return(protocoldata)
}


construct.Illumina = function(sampleSheet=NULL,
			  arrayNames=NULL,
			  ids=NULL,
			  path=".",
			  arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			  highDensity=FALSE,
			  sep="_",
			  fileExt=list(green="Grn.idat",
			  red="Red.idat"),
		      	  cdfName,
		      	  copynumber=TRUE,
		      	  verbose=TRUE, batch, fns, saveDate=TRUE){
       if(!is.null(arrayNames)) {
               pd = new("AnnotatedDataFrame", data = data.frame(Sample_ID=arrayNames))
       }
       if(!is.null(sampleSheet)) { # get array info from Illumina's sample sheet
	       if(is.null(arrayNames)) {
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
	       }
	       pd = new("AnnotatedDataFrame", data = sampleSheet)
	       sampleNames(pd) = basename(arrayNames)
       }
       if(is.null(arrayNames)) {
               arrayNames = gsub(paste(sep, fileExt$green, sep=""), "", dir(pattern=fileExt$green, path=path))
               if(!is.null(sampleSheet)) {
                      sampleSheet=NULL
                      cat("Could not find required info in \'sampleSheet\' - ignoring.  Check \'sampleSheet\' and/or \'arrayInfoColNames\'\n")
               }
               pd = new("AnnotatedDataFrame", data = data.frame(Sample_ID=arrayNames))
       }

       narrays = length(arrayNames)
          
       if(!missing(batch)) {
		stopifnot(length(batch) == narrays)
       }
       if(missing(batch)) {
                batch = as.factor(rep(1, narrays))
       }

       grnfiles = paste(arrayNames, fileExt$green, sep=sep)
       redfiles = paste(arrayNames, fileExt$red, sep=sep)
       if(length(grnfiles)==0 || length(redfiles)==0)
	       stop("Cannot find .idat files")
       if(length(grnfiles)!=length(redfiles))
	       stop("Cannot find matching .idat files")
       if(path[1] != "." & path[1] != ""){
	       grnidats = file.path(path, grnfiles)
	       redidats = file.path(path, redfiles)
       }  else {
	       message("path arg not set.  Assuming files are in local directory, or that complete path is provided")
	       grnidats = grnfiles
	       redidats = redfiles
       }
       if(!all(file.exists(grnidats))) stop("Missing some of the *Grn.idat files")
       if(!all(file.exists(redidats))) stop("Missing some of the *Red.idat files")

	if(verbose) message("Initializing container for genotyping and copy number estimation")
	featureData = getFeatureData.Affy(cdfName, copynumber=copynumber)
	if(!missing(fns)){
		index = match(fns, featureNames(featureData))
		if(all(is.na(index))) stop("fns not in featureNames")
		featureData = featureData[index, ]
	}
	nr = nrow(featureData); nc = narrays
	cnSet = new("CNSet",
		     alleleA=initializeBigMatrix(name="A", nr, nc),
		     alleleB=initializeBigMatrix(name="B", nr, nc),
		     call=initializeBigMatrix(name="call", nr, nc),
		     callProbability=initializeBigMatrix(name="callPr", nr,nc),
		     annotation=cdfName,
		     batch=batch)
        if (!is.null(sampleSheet) && !is.null(sampleSheet$Sample_ID)) {
	   	sampleNames(cnSet) = sampleSheet$Sample_ID
        }
        else sampleNames(cnSet) = arrayNames

	if(saveDate){
		protocolData = getProtocolData.Illumina(grnidats, sep=sep, fileExt=fileExt$green)
	} else{
		protocolData = annotatedDataFrameFrom(A(cnSet), byrow=FALSE)
	}
	rownames(pData(protocolData)) = sampleNames(cnSet)
	protocolData(cnSet) = protocolData
	featureData(cnSet) = featureData
	featureNames(cnSet) = featureNames(featureData)
	pd = data.frame(matrix(NA, nc, 3), row.names=sampleNames(cnSet)) 
	colnames(pd)=c("SKW", "SNR", "gender")
	phenoData(cnSet) = new("AnnotatedDataFrame", data=pd)
	return(cnSet)
}


genotype.Illumina = function(sampleSheet=NULL,
			  arrayNames=NULL,
			  ids=NULL,
			  path=".",
			  arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			  highDensity=FALSE,
			  sep="_",
			  fileExt=list(green="Grn.idat",
			  red="Red.idat"),
		      	  cdfName,
		      	  copynumber=TRUE,
                          batch,
                          fns,
                          saveDate=TRUE,
       			  stripNorm=TRUE,
			  useTarget=TRUE,
		          mixtureSampleSize=10^5,
                          fitMixture=TRUE,
		          eps=0.1,
		          verbose=TRUE,
		          seed=1,
		          sns,
		          probs=rep(1/3, 3),
		          DF=6,
		          SNRMin=5,
		          recallMin=10,
		          recallRegMin=1000,
		          gender=NULL,
		          returnParams=TRUE,
		          badSNP=0.7) {
	is.lds = ifelse(isPackageLoaded("ff"), TRUE, FALSE)
	if(missing(cdfName)) stop("must specify cdfName")
	if(!isValidCdfName(cdfName)) stop("cdfName not valid.  see validCdfNames")
        pkgname = getCrlmmAnnotationName(cdfName)
	callSet = construct.Illumina(sampleSheet=sampleSheet, arrayNames=arrayNames,
			     ids=ids, path=path, arrayInfoColNames=arrayInfoColNames,
                             highDensity=highDensity, sep=sep, fileExt=fileExt, 
			     cdfName=cdfName, copynumber=copynumber, verbose=verbose, batch=batch,
                             fns=fns, saveDate=saveDate)
        if(missing(sns)) sns = sampleNames(callSet)
	open(callSet)
	is.snp = isSnp(callSet)
	snp.index = which(is.snp)
        narrays = ncol(callSet)
        
        if(is.lds) {
          sampleBatches = splitIndicesByNode(seq(along=sampleNames(callSet)))

          mixtureParams = initializeBigMatrix("crlmmMixt-", 4, narrays, "double")
          SNR = initializeBigVector("crlmmSNR-", narrays, "double")
          SKW = initializeBigVector("crlmmSKW-", narrays, "double")
          
          ocLapply(sampleBatches, processIDAT, sampleSheet=sampleSheet, arrayNames=arrayNames,
                 ids=ids, path=path, arrayInfoColNames=arrayInfoColNames, highDensity=highDensity,
                 sep=sep, fileExt=fileExt, saveDate=saveDate, verbose=verbose, mixtureSampleSize=mixtureSampleSize,
                 fitMixture=fitMixture, eps=eps, seed=seed, cdfName=cdfName, sns=sns, stripNorm=stripNorm,
                 useTarget=useTarget, A=A(callSet), B=B(callSet), SKW=SKW, SNR=SNR,
                 mixtureParams=mixtureParams, is.snp=is.snp, neededPkgs=c("crlmm", pkgname))

          open(SKW)
          open(SNR)
          pData(callSet)$SKW = SKW
          pData(callSet)$SNR = SNR
          close(SNR)
          close(SKW)
        } else {
          mixtureParams = matrix(NA, 4, nrow(callSet))

          RG = readIdatFiles(sampleSheet=sampleSheet, arrayNames=arrayNames,
                       ids=ids, path=path, arrayInfoColNames=arrayInfoColNames,
                       highDensity=highDensity, sep=sep, fileExt=fileExt, saveDate=saveDate)

          XY = RGtoXY(RG, chipType=cdfName)
          rm(RG); gc()
        
          res = preprocessInfinium2(XY, mixtureSampleSize=mixtureSampleSize, fitMixture=TRUE, verbose=verbose,
                               seed=seed, eps=eps, cdfName=cdfName, sns=sns, stripNorm=stripNorm, useTarget=useTarget)
          rm(XY); gc()
          if(verbose) message("Finished preprocessing.")
          np.index = which(!is.snp)        
          A(callSet)[snp.index, ] = res[["A"]]
          B(callSet)[snp.index, ] = res[["B"]]
          if(length(np.index)>0) {
            for (j in 1:ncol(callSet)) {
        	A(callSet)[np.index, ] = res[["cnAB"]]$A
        	B(callSet)[np.index, ] = res[["cnAB"]]$B
        	}
          }
	  SKW = pData(callSet)$SKW = res[["SKW"]]
	  SNR = pData(callSet)$SNR = res[["SNR"]]
	  mixtureParams = res[["mixtureParams"]]
          rm(res)
        }

	FUN = ifelse(is.lds, "crlmmGT2", "crlmmGT")
	## genotyping
	crlmmGTfxn = function(FUN,...){
		switch(FUN,
		       crlmmGT2=crlmmGT2(...),
		       crlmmGT=crlmmGT(...))
	}

        if(is.lds) {
          open(A(callSet))
          open(B(callSet))
          tmpA = initializeBigMatrix(name="A", length(snp.index), narrays)
          tmpB = initializeBigMatrix(name="B", length(snp.index), narrays)
          bb = ocProbesets()*length(sns)*8
	  ffrowapply(tmpA[i1:i2, ] <- A(callSet)[snp.index,][i1:i2, ], X=A(callSet)[snp.index,], BATCHBYTES=bb)
	  ffrowapply(tmpB[i1:i2, ] <- B(callSet)[snp.index,][i1:i2, ], X=B(callSet)[snp.index,], BATCHBYTES=bb)
          close(A(callSet))
          close(B(callSet))
          close(tmpA)
          close(tmpB)
        } else {
          tmpA = A(callSet)[snp.index,]
          tmpB = B(callSet)[snp.index,]
        }
        
	tmp = crlmmGTfxn(FUN,
			  A=tmpA,
			  B=tmpB,
			  SNR=SNR,
			  mixtureParams=mixtureParams,
			  cdfName=cdfName,
			  row.names=NULL,
			  col.names=sampleNames(callSet),
			  probs=probs,
			  DF=DF,
			  SNRMin=SNRMin,
			  recallMin=recallMin,
			  recallRegMin=recallRegMin,
			  gender=gender,
			  verbose=verbose,
			  returnParams=returnParams,
			  badSNP=badSNP)
	if(verbose) message("Genotyping finished.  Updating container with genotype calls and confidence scores.")
	if(is.lds){
		bb = ocProbesets()*ncol(callSet)*8
		open(tmp[["calls"]])
		open(tmp[["confs"]])
		ffrowapply(snpCall(callSet)[snp.index,][i1:i2, ] <- tmp[["calls"]][i1:i2, ], X=tmp[["calls"]], BATCHBYTES=bb)
		ffrowapply(snpCallProbability(callSet)[snp.index,][i1:i2, ] <- tmp[["confs"]][i1:i2, ], X=tmp[["confs"]], BATCHBYTES=bb)
#		close(tmp[["calls"]])
#		close(tmp[["confs"]])
#                open(tmpA); open(tmpB)
#                delete(tmpA); delete(tmpB);
                delete(tmp[["calls"]]); delete(tmp[["confs"]])
                rm(tmpA, tmpB)
	} else {
		calls(callSet)[snp.index, ] = tmp[["calls"]]
		snpCallProbability(callSet)[snp.index, ] = tmp[["confs"]]
                rm(tmpA, tmpB)
	}
	callSet$gender = tmp$gender
        rm(tmp)
	close(callSet)
	return(callSet)
}


processIDAT =  function(sel, sampleSheet=NULL,
			  arrayNames=NULL,
			  ids=NULL,
			  path=".",
			  arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			  highDensity=FALSE,
			  sep="_",
			  fileExt=list(green="Grn.idat", red="Red.idat"),
			  saveDate=FALSE, 
                          verbose=TRUE,
                          mixtureSampleSize=10^5,
			  fitMixture=TRUE,
			  eps=0.1,
			  seed=1,
			  cdfName,
			  sns,
			  stripNorm=TRUE,
			  useTarget=TRUE,
                          A, B, SKW, SNR, mixtureParams, is.snp) {
  
        if(length(path)>= length(sel)) path = path[sel]
        RG = readIdatFiles(sampleSheet=sampleSheet[sel,], arrayNames=arrayNames[sel],
                       ids=ids, path=path, arrayInfoColNames=arrayInfoColNames,
                       highDensity=highDensity, sep=sep, fileExt=fileExt, saveDate=saveDate)

        XY = RGtoXY(RG, chipType=cdfName)
        open(RG@assayData$R); open(RG@assayData$G); open(RG@assayData$zero)
        delete(RG@assayData$R); delete(RG@assayData$G); delete(RG@assayData$zero); rm(RG)
        gc()
        if (missing(sns)) sns = sampleNames(XY)
        
        res = preprocessInfinium2(XY, mixtureSampleSize=mixtureSampleSize, fitMixture=TRUE, verbose=verbose,
                               seed=seed, eps=eps, cdfName=cdfName, sns=sns, stripNorm=stripNorm, useTarget=useTarget)
        #                       save.it=save.it, snpFile=snpFile, cnFile=cnFile)
        open(XY@assayData$X); open(XY@assayData$Y); open(XY@assayData$zero)
        delete(XY@assayData$X); delete(XY@assayData$Y); delete(XY@assayData$zero); rm(XY)
        gc()
	if(verbose) message("Finished preprocessing.")
        snp.index = which(is.snp)
	np.index = which(!is.snp)      

        open(res[["A"]])
        open(res[["B"]])
        open(res[["SKW"]])
        open(res[["SNR"]])
        open(res[["mixtureParams"]])
	bb = ocProbesets()*length(sns)*8
        ffrowapply(A[snp.index,][i1:i2, sel] <- res[["A"]][i1:i2, ], X=res[["A"]], BATCHBYTES=bb)
	ffrowapply(B[snp.index,][i1:i2, sel] <- res[["B"]][i1:i2, ], X=res[["B"]], BATCHBYTES=bb)
	if(length(np.index)>0) {
          for (j in sel) {
            A[np.index, j] = res[["cnAB"]]$A[,j]
            B[np.index, j] = res[["cnAB"]]$B[,j]
          }
        }
        delete(res[["A"]]); delete(res[["B"]])
	SKW[sel] = res[["SKW"]][]
	SNR[sel] = res[["SNR"]][]
	mixtureParams[,sel] = res[["mixtureParams"]][]
        close(A)
        close(B)
        close(SNR)
        close(SKW)
        close(mixtureParams)
        delete(res[["SKW"]]); delete(res[["SNR"]]); delete(res[["mixtureParams"]])
        rm(res)
        TRUE
      }
