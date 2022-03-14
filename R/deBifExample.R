#' Examples of phaseplane analysis of a system of ODEs
#'
#' \code{deBifExample}
#'
#'
#'   deBifExample(example)
#'
#'
#' Function runs one of the examples provided with the deBif package
#'
#' @param   example  (string, optional)
#' \preformatted{}
#'               Name of the example. If not provided a list of examples
#'               is returned
#'
#' @return None.
#'
#' @export
deBifExample <- function(example) {
  # locate all the shiny app examples that exist
  validExamples <- gsub(".R$", "", list.files(system.file("shiny-examples", package = "deBif")))

  validExamplesMsg <-
    paste0(
      "Valid examples are: '",
      paste(validExamples, collapse = "', '"),
      "'")

  if (!missing(example)) example <- gsub(".R$", "", example)

  # if no example name is given, throw an error
  if (missing(example) || !nzchar(example) || (!(example %in% validExamples))) {
    stop(
      'Please run `deBifExample()` with a valid example app as an argument.\n',
      validExamplesMsg,
      call. = FALSE)
  }

  # find and launch the app
  appDir <- system.file("shiny-examples", paste0(example, ".R"), package = "deBif")
  # file.show(appDir)
  shiny::runApp(appDir, display.mode = "normal")
}
