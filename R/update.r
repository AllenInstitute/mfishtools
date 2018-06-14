#' Update the mfishtools library
#' 
update_mfishtools <- function() {
  devtools::install_github("AllenInstitute/mfishtools",
    auth_token="802976690281f1483c40de46d0a07e9d01a3de08")
}