#' Check that the right SQLite database was given
#'
#' This function checks the suufix of the SQLite database path provided in the
#' outer function and halts if it does not match the name on MethMotif server.
#' This is an internal helper function not meant to be called on its own.
#' @param local_db_path TFregulome url provided to the outer function
#' @return
#' @keywords internal
#' @export
#' @examples
#' check_db_file(local_db_path)
#'

check_db_file <- function(local_db_path = NULL) {
  # check file name
  if (!is.null(local_db_path)) {
    if (endsWith(local_db_path, suffix = "/tfregulome.sqlite")==FALSE) {
      stop("local SQLite database should be call 'tfregulome.sqlite'!")
    } else if (!file.exists(local_db_path)) {
      stop("cannot locate the SQLite database from the provided path!")
    }
  }
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

query_local_database <- function(db, query_index, query_value, id) {
  # get the root path of the local tfregulome directory
  base_dir <- paste0(dirname(db), "/TFregulome_database")

  # prepare query
  if (!missing(id)) {
    query <- paste0("SELECT * FROM TFBS_table WHERE UPPER(id)=UPPER('", id, "')")
  }
  else if (!missing(query_index) && !missing(query_value)) {
    if (sum(query_index) == 0) {
      query <- "SELECT * FROM TFBS_table"
    }
    else {
      query <- paste0("SELECT * FROM TFBS_table WHERE ",
                      paste0("UPPER(", sub("=", ")=UPPER('", query_value[query_index == 1]),
                             "')", collapse = " AND "))
    }
  }
  else {
    stop("Trying to query database without parameters")
  }

  # query the database
  con <- dbConnect(RSQLite::SQLite(), db)
  results <- tryCatch({
    dbGetQuery(con, query)
  },
  error = function(cond) {
    message("There is a warning to connect TFregulomeR SQLite database!")
    message("Advice:")
    message("1) Check the path to the local database;")
    message("2) Check dependent package 'RSQLite';")
    message(paste0("warning: ",cond))
    return(NULL)
  })
  dbDisconnect(con)

  if (nrow(results)>0) {
    # modify the output table
    results$motif_MEME <- paste0(base_dir, "_", results$species, "/",
                                 results$organ, "/motif_matrix/",
                                 results$motif_MEME)
    results$motif_TRANSFAC <- paste0(base_dir, "_", results$species, "/",
                                     results$organ, "/motif_matrix/",
                                     results$motif_TRANSFAC)
    results$beta_score_matrix <- paste0(base_dir, "_", results$species, "/",
                                        results$organ, "/beta_score_matrix/",
                                        results$beta_score_matrix)
    results$all_peak_file <- paste0(base_dir, "_", results$species, "/",
                                    results$organ, "/TF_all_peaks/",
                                    results$all_peak_file)
    results$peak_with_motif_file <- paste0(base_dir, "_", results$species, "/",
                                           results$organ, "/TF_peaks_with_motif/",
                                           results$peak_with_motif_file)
    results$DNA_methylation_profile <- paste0(base_dir, "_", results$species, "/",
                                              results$organ, "/DNA_methylation_profile/",
                                              results$DNA_methylation_profile)
    results$DNA_methylation_profile_200bp <- paste0(base_dir, "_",
                                                    results$species, "/",
                                                    results$organ,
                                                    "/DNA_methylation_profile_200bp/",
                                                    results$DNA_methylation_profile_200bp)
    results$TFBS <- paste0(base_dir, "_", results$species, "/", results$organ,
                           "/TFBS/", results$TFBS)
    results$Ncor_between_MEME_ChIP_and_HOMER <- as.logical(results$Ncor_between_MEME_ChIP_and_HOMER)

    return(results)
  }
  else {
    return(results)
  }
}
