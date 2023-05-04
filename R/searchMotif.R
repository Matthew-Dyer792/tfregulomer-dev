#' Search motif PFM and beta score matrix (if source is MethMotif) for a given TFregulomeR ID in TFregulomeR
#'
#' This function allows you to obtain motif PFM matrix and beta score matrix (if source is MethMotif) for a given TFregulomeR ID in TFregulomeR
#' @param id Required. TFregulomeR ID.
#' @param motif_format Motif PFM format, either in MEME by default or TRANSFAC.
#' @param server server localtion to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "https://bioinfo-csi.nus.edu.sg/methmotif/" or "https://methmotif.org", please use a new url.
#' @param local_db_path The complete path to the SQLite implementation of TFregulomeR database available at "https://methmotif.org/API_ZIPPED.zip"
#' @return MethMotif class object
#' @keywords MethMotif
#' @export
#' @examples
#' K562_CEBPB <- searchMotif(id = "MM1_HSA_K562_CEBPB")
#' K562_CEBPB_transfac <- searchMotif(id = "MM1_HSA_K562_CEBPB",
#'                                    motif_format = "TRANSFAC")

searchMotif <- function(id, motif_format = "MEME",
                        server = "ca", TFregulome_url, local_db_path = NULL)
{
  # check motif_format MEME and TRANSFAC.
  motif_format = toupper(motif_format)
  if (motif_format != "MEME" & motif_format != "TRANSFAC")
  {
    stop("Please check motif_format! Currently we only support MEME (default) and TRANSFAC formats!")
  }

  # call API helper function
  TFregulome_url <- construct_API_url(server, TFregulome_url)
  # helper function to check SQLite database
  check_db_file(local_db_path)

  if(missing(id))
  {
    stop("Please input regulome ID using 'id = '!")
  }
  else if (!is.null(local_db_path)) {
    # make a request to the local database
    request_content_df <- query_local_database(local_db_path, id = id)
  }
  else
  {
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
    message("There are a matched record exported in a MethMotif object.")
    # using the previous backup TFregulome_url to parse into readMMmotif()
    if (motif_format == "MEME")
    {
      motif_file_path <- request_content_df[1, "motif_MEME"]
    }
    else
    {
      motif_file_path <- request_content_df[1, "motif_TRANSFAC"]
    }
    betaScore_file_path <- request_content_df[1, "beta_score_matrix"]
    num_peak <- request_content_df[1,"peak_with_motif_num"]
    nsites <- request_content_df[1,"TFBS_num"]

    # MEME file path for background extraction when motif formati is TRANSFAC
    motif_file_path_MEME <- request_content_df[1, "motif_MEME"]
    methmotif_item_motif <- readMMmotif(motif_file_path = motif_file_path,
                                        motif_format = motif_format,
                                        id = id,
                                        num_peak = as.integer(num_peak),
                                        nsites = as.integer(nsites),
                                        motif_file_path_MEME = motif_file_path_MEME)
    methmotif_item_betaScore <- readBetaScoreMatrix(betaScore_file_path = betaScore_file_path)
    methmotif_item <- new("MethMotif")
    methmotif_item <- updateMethMotif(methmotif_item,
                                      MMBetaScore = methmotif_item_betaScore,
                                      MMmotif = methmotif_item_motif)

    return(methmotif_item)
  }
}
