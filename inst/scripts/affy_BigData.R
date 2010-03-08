## We need to initialize the data storage mode outside of the core
library(crlmm)
library(ff)
outdir <- "/amber1/scratch/rscharpf/hapmap/affy_bigdata"
options(fffinalizer="close", fftempdir=outdir, fffinonexit="TRUE")
PATH <- "/thumper/ctsa/snpmicroarray/hapmap/raw/affy/phase3_1m"
celFiles <- list.celfiles(PATH, full.names=TRUE, pattern=".CEL")
obj <- constructClass("genomewidesnp6")
crlmmOptions(obj)
##
##
##	obj <- construct(obj, filenames)
##	obj <- snprma(obj, filenames)
##	ops <- crlmmOptions(obj)
##	if(ops$copynumber & ops$nonpolymorphic)
##		obj <- cnrma(obj, filenames)
##	message("Initializing CallSet")
##	object <- as(obj, "CallSet")
##	rm(obj); gc()
##	callSet <- crlmm.batch(object)
##
callSet <- crlmm(obj, celFiles)
q("no")
##save(callSet, file=file.path(outdir, "callSet.rda"))
##load(file.path(outdir, "callSet.rda"))
callSet$batch <- substr(basename(celFiles), 13, 13)
## Default options for CN estimation
crlmmOptions(callSet)$cnOpts
tmp <- as(callSet, "CNSet")
tmp2 <- tmp[1:100, 1:10]
cnSet <- trace(computeCopynumber.CNSet, browser)




