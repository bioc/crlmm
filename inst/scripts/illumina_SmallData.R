library(crlmm)
outdir <- "/amber1/scratch/rscharpf/hapmap/illumina_smalldata"
obj <- constructClass("human370v1c")
##outdir <- "/amber1/scratch/rscharpf/test/illumina_SmallData"
datadir <- "/thumper/ctsa/snpmicroarray/illumina/IDATS/370k"
samplesheet = read.csv(file.path(datadir, "HumanHap370Duo_Sample_Map.csv"), header=TRUE, as.is=TRUE)
samplesheet <- samplesheet[-c(28:46,61:75,78:79), ]
filenames <- file.path(datadir, unique(samplesheet[, "SentrixPosition"]))
crlmmOptions(obj)$readOpts$sampleSheet <- samplesheet
callSet <- crlmm(obj, filenames)



	RG <- construct(obj, filenames)
	RG <- readIdatFiles(RG, filenames)
	XY <- RGtoXY(RG)
	rm(RG); gc()
	storageMode(XY) <- "environment"
	stripNorm <- crlmmOptions(XY)$readOpts[["stripNorm"]]
	useTarget <- crlmmOptions(XY)$readOpts[["useTarget"]]
	verbose <- crlmmOptions(XY)$verbose
	if(stripNorm)
		XY = stripNormalize(XY, useTarget=useTarget, verbose=verbose)
	alleleSet <- as(XY, "IlluminaAlleleSet")
	##Do i need to pass zero.  Is the only thing updated the mixtureParams and SNR?
	alleleSet <- preprocessInfinium2(alleleSet, zero=Z(XY))

	callSet <- as(alleleSet, "CallSet")
trace(crlmm.CallSet, browser)
crlmm.batch(callSet)

