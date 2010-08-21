setOldClass("ellipse")
setOldClass("ff_matrix")
setOldClass("ffdf")
setClassUnion("ff_or_matrix", c("ff_matrix", "matrix", "ffdf"))

