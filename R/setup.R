.onAttach <- function(libname, pkgname) {

  msg <- "\nWelcome to the deBif package for analysis of ordinary differential equations\nExecute the deBifHelp() function for instructions\nUse deBifReset() to reload the package in case of computational problems\n"

  options(stringsAsFactors=FALSE)
  packageStartupMessage(msg)
}

.onUnload <- function(libpath) {
  library.dynam.unload("deBif", libpath)
}
