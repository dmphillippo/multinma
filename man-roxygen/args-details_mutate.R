#' @details All arguments specifying columns of `data` accept the following:
#'
#'  * A column name as a character string, e.g. `study = "studyc"`
#'  * A bare column name, e.g. `study = studyc`
#'  * `dplyr::mutate()` style semantics for inline variable transformations, e.g. `study = paste(author, year)`
#'
