#' convert motif PFM in TFregulomeR into PFMatrix class object in TFBSTools package
#'
#' This function allows you to retrieve and convert motif PFM in TFregulomeR database into PFMatrix class object, which can be used in TFBSTools package.
#' @param id Required. TFregulomeR ID.
#' @param server server localtion to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "https://bioinfo-csi.nus.edu.sg/methmotif/" or "https://methmotif.org", please use a new url.
#' @param local_db_path The complete path to the SQLite implementation of TFregulomeR database available at "https://methmotif.org/API_ZIPPED.zip"
#' @return  An object of class PFMatrix
#' @keywords toTFBSTools
#' @export
#' @examples
#' require(TFBSTools)
#' CEBPB_pfm <- toTFBSTools(id = "MM1_HSA_K562_CEBPB")

toTFBSTools <- function(id, server = "ca", TFregulome_url,
                        local_db_path = NULL) {
  if (missing(id)) {
    stop("Please provide a TFregulomeR ID using 'id ='")
  }

  # build api_object
  api_object <- .construct_api(server, TFregulome_url, local_db_path)
  methmotif_output <- suppressMessages(
    .searchMotif(
      id = id,
      motif_format = "TRANSFAC",
      api_object = api_object
    )
  )
  if (is.null(methmotif_output)) {
    message(paste0("No record for id ", id, " in TFregulomeR!"))
    return(NULL)
  } else {
    methmotif_output_transfac <- methmotif_output
    pfm <- TFBSTools::PFMatrix(
      ID = methmotif_output_transfac@MMmotif@id,
      name = methmotif_output_transfac@MMmotif@alternate_name, strand = "*",
      bg = methmotif_output_transfac@MMmotif@background,
      profileMatrix = t(methmotif_output_transfac@MMmotif@motif_matrix)
    )
    return(pfm)
  }
}
