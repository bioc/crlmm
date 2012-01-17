setOldClass("ellipse")
setOldClass("ffdf")
##setClassUnion("ff_or_matrix", c("ffdf", "ff_matrix", "matrix"))
setClass("PredictionRegion", contains="list")
