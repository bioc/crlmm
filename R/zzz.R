# Loading required libraries
THISPKG <- "crlmm"

.onLoad <- function(libname, pkgname) {
  require("methods")
}

.onAttach <- function(libname, pkgname) {
  message("Welcome to crlmm version ", packageDescription(THISPKG, field="Version"))
}

.onUnload <- function( libpath ){
  library.dynam.unload(THISPKG, libpath)
}

  .crlmmPkgEnv <- new.env(parent=emptyenv())
