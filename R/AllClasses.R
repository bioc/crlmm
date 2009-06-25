## Class definition
setClass("ABset", contains="eSet",
	 prototype = prototype(
	 new("VersionedBiobase",
	     versions=c(classVersion("eSet"), SnpSet="1.0.0"))))
setClass("crlmmSet", contains="eSet")
setClass("CrlmmSetList", contains="list")
setClass("CopyNumberSet", contains="eSet",
	 prototype = prototype(
	 new("VersionedBiobase",
	     versions=c(classVersion("eSet"), SnpSet="1.0.0"))))
