#' Standardized api object constructor
#'
#' This function creates a standardized function to check if the TFregulome
#' server is valid. If that is successful it will construct a api-object and
#' return it. This is an internal helper function not meant to be called on
#' its own.
#' @param server server location from the calling function
#' @param TFregulome_url specified url from the calling function
#' @param local_db_path path to the local db from the calling function
#' @return api_object
#' @keywords internal
#' @examples
#' .construct_api(server = "ca")
#' .construct_api(TFregulome_url = "https://methmotif.org/")
#' .construct_api(local_db_path = "/opt/database/tfregulome.sqlite")
#'

.construct_api <- function(server, TFregulome_url, local_db_path) {
  # build api_object
  if (!is.null(local_db_path)) {
    api_object <- LocalDatabase(address = local_db_path)
  } else if (server != "sg" && server != "ca") {
    stop("server should be either 'ca' (default) or 'sg'!")
  } else {
    api_object <- WebDatabase(address = TFregulome_url, server = server)
  }
  return(api_object)
}
