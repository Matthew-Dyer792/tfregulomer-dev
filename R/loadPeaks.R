#' load peaks from TFregulomeR
#'
#' This function allows you to obtain the peaks from TFregulomeR using TFregulomeR ID.
#' @param id Required. TFregulomeR ID
#' @param includeMotifOnly Either TRUE or FALSE (default). If TRUE, only peaks with motif will be returned
#' @param server server localtion to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "https://bioinfo-csi.nus.edu.sg/methmotif/" or "https://methmotif.org", please use a new url.
#' @param local_db_path The complete path to the SQLite implementation of TFregulomeR database available at "https://methmotif.org/API_TFregulomeR/downloads/"
#' @return  a data.frame containing peak coordinates
#' @keywords loadPeaks
#' @export
#' @examples
#' CEBPB_peaks <- loadPeaks(id = "MM1_HSA_K562_CEBPB")

loadPeaks <- function(id, includeMotifOnly = FALSE,
                      server = "ca", TFregulome_url, local_db_path = NULL)
{
  # call API helper function
  TFregulome_url <- construct_API_url(server, TFregulome_url)
  # helper function to check SQLite database
  check_db_file(local_db_path)

  # check input argument id
  if (missing(id)) {
    stop("Please input a TFregulomeR id using 'id = '.")
  }
  else if (!is.null(local_db_path)) {
    # make a request to the local database
    request_content_df <- query_local_database(local_db_path, id = id)
  }
  else {
    # make a json request to the API
    request_content_json <- API_request(TFregulome_url, id = id)
    if (is.null(request_content_json)) {
      message("Empty output for the request!")
      return(NULL)
    }
    else {
      request_content_df <- as.data.frame(request_content_json$TFBS_records)
    }
  }
  # process the returned data.frame
  if (nrow(request_content_df)==0) {
    if (exists("request_content_json")) {
      message(request_content_json$message)
      return(NULL)
    } else {
      message("No matched records found.")
      return(NULL)
    }
  }
  else {
    if (includeMotifOnly) {
      peak_file <- request_content_df[1,c("peak_with_motif_file")]
    }
    else {
      peak_file <- request_content_df[1,c("all_peak_file")]
    }
    # read peak file
    peak_df <- tryCatch({
      read.delim(peak_file, sep = "\t", header = FALSE)
    },
    warning = function(w) {
      message(paste0("Warning: No peak file for id=", id))
      stop()
    },
    error = function(err) {
      message(paste0("Error: No peak file for id=", id))
    })
    colnames(peak_df) <- c("chr","start","end","id","tag_fold_change")
    # use new id for each peak region
    if (includeMotifOnly) {
      peak_df$id = paste0(id,"_peaks_with_motif_", as.vector(rownames(peak_df)))
    }
    else {
      peak_df$id = paste0(id,"_all_peaks_", as.vector(rownames(peak_df)))
    }
    message("Success: peak file has been returned in a data frame!")
    return(peak_df)
  }
}
