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
                      server = "ca", TFregulome_url,
                      local_db_path = NULL) {
  # build api_object
  api_object <- .construct_api(server, TFregulome_url, local_db_path)
  results <- .loadPeaks(
    id = id,
    includeMotifOnly = includeMotifOnly,
    api_object = api_object
  )
  return(results)
}

#' load peaks from TFregulomeR
#'
#' This function allows you to obtain the peaks from TFregulomeR using TFregulomeR ID.
#' @param id Required. TFregulomeR ID
#' @param includeMotifOnly Either TRUE or FALSE (default). If TRUE, only peaks with motif will be returned
#' @param api_object s4 object for interacting with the TFregulome database
#' @return  a data.frame containing peak coordinates
#' @keywords loadPeaks
#' @examples
#' CEBPB_peaks <- .loadPeaks(id = "MM1_HSA_K562_CEBPB")

.loadPeaks <- function(id, includeMotifOnly = FALSE, api_object) {
  # check input argument id
  if (missing(id)) {
    stop("Please input a TFregulomeR id using 'id = '.")
  } else if (!is(api_object, "API") && !validObject(api_object)) {
    stop("Invalid API object!")
  } else {
    # make the request
    request_content_df <- apiRequest(api_object, id = id)
  }
  # check db output
  if (is.null(request_content_df)) {
    return(NULL)
  } else {
    if (includeMotifOnly) {
      peak_file <- request_content_df[1, c("peak_with_motif_file")]
    } else {
      peak_file <- request_content_df[1, c("all_peak_file")]
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
    colnames(peak_df) <- c("chr", "start", "end", "id", "tag_fold_change")
    # use new id for each peak region
    if (includeMotifOnly) {
      peak_df$id <- paste0(
        id, "_peaks_with_motif_",
        as.vector(rownames(peak_df))
      )
    } else {
      peak_df$id <- paste0(id, "_all_peaks_", as.vector(rownames(peak_df)))
    }
    message("Success: peak file has been returned in a data frame!")
    return(peak_df)
  }
}