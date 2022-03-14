.onAttach <- function(libname, pkgname) {

  msg <- "\nWelcome to the deBif package for analysis of ordinary differential equations\nExecute the deBifHelp() function for instructions\n"

  options(stringsAsFactors=FALSE)
  packageStartupMessage(msg)
}
