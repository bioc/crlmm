library(crlmm)
library(ff)
outdir <- "/amber1/scratch/rscharpf/hapmap/illumina_bigdata"
if(!file.exists(outdir)) dir.create(outdir)
options(fffinalizer="close", fftempdir=outdir, fffinonexit="TRUE")
obj <- constructClass("human370v1c")
##outdir <- "/amber1/scratch/rscharpf/test/illumina_SmallData"
datadir <- "/thumper/ctsa/snpmicroarray/illumina/IDATS/370k"
samplesheet = read.csv(file.path(datadir, "HumanHap370Duo_Sample_Map.csv"), header=TRUE, as.is=TRUE)
samplesheet <- samplesheet[-c(28:46,61:75,78:79), ]
filenames <- file.path(datadir, unique(samplesheet[, "SentrixPosition"]))
crlmmOptions(obj)$readOpts$sampleSheet <- samplesheet
callSet <- crlmm(obj, filenames)



