setClass("CopyNumberSet", contains="SnpLevelSet")
setClass("CrlmmSet", contains=c("SnpCallSetPlus", "CopyNumberSet"))
setClass("SnpCallSetPlusFF", contains="SnpCallSetPlus")
setClass("CrlmmSetFF", contains="CrlmmSet")
setClass("SegmentSet", contains="CrlmmSet",
	 representation(emissionPr="array",
			segmentData="data.frame"))
