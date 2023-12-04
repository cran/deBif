#' Reloads the deBif package
#'
#' \code{deBifReset} unloads and reloads the deBif package.
#'
#' @return None.
#'
#' @examples
#' if(interactive()){
#' deBifReset()
#' }
#'
#' @export
deBifReset <- function() {
  unloadNamespace("deBif")
  suppressMessages(require(deBif, quietly = TRUE))
}

