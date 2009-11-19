##setClass("CNSet", contains="eSet")
##setClass("CNSet", contains=c("SnpCallSetPlus", "CNSet"))
setClass("CNSet", contains="SnpCallSetPlus",
	 representation(emissionPr="array",
			segmentData="RangedData"))

##setClass("SegmentSet", contains="CNSet",
##	 representation(emissionPr="array",
##			segmentData="data.frame"))
