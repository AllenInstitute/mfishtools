#' Update the mfishtools library
#'
#' @export
update_mfishtools <- function() {
  devtools::install_github("AllenInstitute/mfishtools", build_vignettes = TRUE)
}
