krlmm = function(cnSet, cdfName, gender) { # Doesn't do anything as yet - Jenny/Cynthia to add functionality based on their code
  if(is.null(gender)) {
    gender=rep(1, ncol(cnSet)) # Need to impute gender if not specified - Cynthia is working on this using Y chr SNPs
  }
  open(cnSet$gender)
  cnSet$gender[,] = gender
  close(cnSet$gender)

  TRUE
}
