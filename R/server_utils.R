#' Standardized format for API url
#'
#' This function creates a standardized TFregulome url in order to prevent
#' connection errors. This is an internal helper function not meant to be
#' called on its own.
#' @param server location provided to the outer function
#' @param TFregulome_url TFregulome url provided to the outer function
#' @return properly formatted TFregulome url
#' @keywords internal
#' @export
#' @examples
#' construct_API_url(server, TFregulome_url)
#'

construct_API_url <- function(server = "ca", TFregulome_url) {
  # check server location
  if (server != "sg" && server != "ca") {
    stop("server should be either 'ca' (default) or 'sg'!")
  }

  # make an appropriate API url
  if (missing(TFregulome_url)){
    if(server == 'sg') {
      TFregulome_url <- "https://bioinfo-csi.nus.edu.sg/methmotif/api/table_query/"
    } else {
      TFregulome_url <- "https://methmotif.org/api/table_query/"
    }
  } else if (endsWith(TFregulome_url, suffix = "/index.php")==TRUE){
    TFregulome_url <- gsub("index.php", "", TFregulome_url)
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else if (endsWith(TFregulome_url, suffix = "/")==TRUE){
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else {
    TFregulome_url <- paste0(TFregulome_url, "/api/table_query/")
  }

  return(TFregulome_url)
}


#' Browse the current data available in TFregulomeR
#'
#' This function allows you to get the current data in TFregulomeR
#' @param db A complete path to the SQLite database
#' @param query The query statement to be run on the SQLite database
#' @return  data.frame containing the information of the queried TFBSs in the local TFregulomeR database
#' @keywords internal SQLite
#' @export
#' @examples
#' request_content_df <- query_local_database("/home/USER/tfregulome.sqlite", "SELECT * FROM TFBS_table")

API_request <- function(TFregulome_url, query_index, query_value, id) {
  # prepare query_url
  if (!missing(id)) {
    query_url <- paste0("listTFBS.php?AllTable=F&id=", id)
  } else if (!missing(query_index) && !missing(query_value)) {
    if (sum(query_index) == 0) {
      query_url <- "listTFBS.php?AllTable=T"
    } else {
      query_url <- paste0("listTFBS.php?AllTable=F&",
                          paste0(query_value[query_index == 1], collapse = "&"))
    }
  } else {
    stop("Trying to query database without parameters")
  }

  #parse JSON from API endpoint
  request_content_json <- tryCatch({
    fromJSON(paste0(TFregulome_url, query_url))
  },
  error = function(cond) {
    message("There is a warning to connect TFregulomeR API!")
    message("Advice:")
    message("1) Check internet access;")
    message("2) Check dependent package 'jsonlite';")
    message("3) Current TFregulomeR server is implemented in MethMotif database, whose homepage is 'https://bioinfo-csi.nus.edu.sg/methmotif/' or 'https://methmotif.org'. If MethMotif homepage url is no more valid, please Google 'MethMotif', and input the valid MethMotif homepage url using 'TFregulome_url = '.")
    message(paste0("warning: ", cond))
    return(NULL)
  })
  return(request_content_json)
}
