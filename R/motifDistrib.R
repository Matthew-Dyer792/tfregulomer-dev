#' profile TFBS distribution
#'
#' This function allows you to profile TFBS distributions in a given list of peak sets.
#' @param id Required. TFregulomeR ID. The TFBS of interest to be profiled.
#' @param peak_list Required. List of data.frames, each of which contains bed-format peak regions. They are the peak sets in which you want to profile the TFBS distributions, and can be loaded from TFregulomeR database or self-provided.
#' @param peak_id Required. Character of vector, each of which is a unique ID corresponding to the element in "peak_list". If a peak set is orignally from TFregulomeR Database, its TFregulomeR ID should be used here.
#' @param plot_at_each_side By default 100bp, and motif occurrences in a window of +/- 100bp around peak centres will be returned.
#' @param server server localtion to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "https://bioinfo-csi.nus.edu.sg/methmotif/" or "https://methmotif.org", please use a new url.
#' @param local_db_path The complete path to the SQLite implementation of TFregulomeR database available at "https://methmotif.org/API_ZIPPED.zip"
#' @return  a list containing the numbers of input peaks and peaks with motif, as well as motif occurrences in the plot window.
#' @keywords motifDistrib
#' @export
#' @examples
#' CEBPB_peaks <- loadPeaks(id = "MM1_HSA_K562_CEBPB")
#' motifDistrib_output <- motifDistrib(id = "MM1_HSA_K562_CEBPB",
#'                                     peak_list = list(CEBPB_peaks[1:100,]),
#'                                     peak_id = "MM1_HSA_K562_CEBPB")

motifDistrib <- function(id, peak_list, peak_id, plot_at_each_side = 100,
                         server = "ca", TFregulome_url, local_db_path = NULL) {
  # check input arguments
  if (missing(id)) {
    stop("Please provide a TFregulomeR ID using 'id ='!")
  }
  if (missing(peak_list)) {
    stop("Please provide peak sets stored in list class list() using 'peak_list ='! Peak sets can be loaded from TRregulomeR using loadPeaks() or your own peaks in bed-like format")
  }
  if (missing(peak_id)) {
    stop("Please provide unique ids in vector corresponding to the peak sets in 'peak_list'. If the peak set was derived from TFregulomeR, please use its TFregulomeR ID here!")
  }
  if (length(peak_list) != length(peak_id)) {
    stop("The number of peak sets in 'peak_list' is NOT equal to the number of ids in 'peak_id'")
  }
  if (!is.character(id)) {
    stop("id should be class character!")
  }
  if (!is.list(peak_list)) {
    stop("'peak_list' should be class list!")
  }
  if (!is.numeric(plot_at_each_side)) {
    stop("'plot_at_each_side' should be class numeric!")
  }

  # build api_object
  api_object <- .construct_api(server, TFregulome_url, local_db_path)

  # start analysing
  message(paste0("motifDistrib starts analysing for TFregulomeR ID = ", id))
  # get motif sequences for id
  request_content_df <- apiRequest(api_object,  id = id)
  # process the returned data.frame
  if (nrow(request_content_df) == 0) {
    message(paste0("No record was found for your input TFregulomeR ID. Your input: id = ", id, "."))
    return(NULL)
  } else {
    motif_seq_path <- request_content_df[1, c("TFBS")]
    motif_seq <- read.delim(motif_seq_path, sep = "\t", header = FALSE)
    colnames(motif_seq) <- c("motif_chr", "motif_start", "motif_end",
                              "motif_strand", "motif_weight", "motif_pvalue",
                              "motif_qvalue", "motif_sequence")
    motif_seq$motif_id <- paste0(
      id, "_motif_sequence_", as.vector(rownames(motif_seq))
    )
    motif_seq_grange <- GRanges(
      motif_seq$motif_chr,
      IRanges(motif_seq$motif_start + 1,
              motif_seq$motif_end),
      names = motif_seq$motif_id
    )
    names(motif_seq_grange) <- motif_seq$motif_id
  }
  motifDistrb_list <- list()
  for (i in seq(1, length(peak_id), 1)) {
    peak_id_i <- peak_id[i]
    message(paste0("... ... analysing peak set ", peak_id_i))
    peak_target <- peak_list[[i]]
    # amend peak set
    peak_target <- peak_target[, c(1, 2, 3)]
    colnames(peak_target) <- c("peak_chr", "peak_start", "peak_end")
    peak_target$target_peak_id <- paste0(
      peak_id_i, "_target_peak_", as.vector(rownames(peak_target))
    )

    # # query MethMotif API
    # if (!is.null(local_db_path)) {
    #   # make a request to the local database
    #   request_content_df <- query_local_database(local_db_path, id = id)
    # } else {
    #   # make a json request to the API
    #   request_content_json <- API_request(TFregulome_url, id = id)
    #   if (is.null(request_content_json)) {
    #     isTFregulomeID <- FALSE
    #   } else {
    #     request_content_df <- as.data.frame(request_content_json$TFBS_records)
    #   }
    # }
    # make the request
    request_content_df <- apiRequest(api_object, id = id)
    # check json content obtained from MethMotif API
    if (nrow(request_content_df) == 0) {
      isTFregulomeID <- FALSE
    } else {
      isTFregulomeID <- TRUE
    }
    if (isTFregulomeID) {
      peak_target_grange <- GRanges(peak_target$peak_chr,
                                    IRanges(peak_target$peak_start - 99,
                                            peak_target$peak_end + 100),
                                    names = peak_target$target_peak_id)
    } else  {
      peak_target_grange <- GRanges(peak_target$peak_chr,
                                    IRanges(peak_target$peak_start,
                                            peak_target$peak_end),
                                    names = peak_target$target_peak_id)
    }
    # overlap with motif seq
    # subsetOverlaps may mis-think the two sets coming from different references, so suppressWarnings here
    suppressWarnings(overlapped <- subsetByOverlaps(peak_target_grange, motif_seq_grange))
    # combine overlapped motif seq information with peak
    hits <- findOverlaps(peak_target_grange, motif_seq_grange)
    motif_id <- CharacterList(split(names(motif_seq_grange)[subjectHits(hits)],
                                             queryHits(hits)))
    mcols(overlapped) <- DataFrame(mcols(overlapped), motif_id)
    overlapped_df <- as.data.frame(overlapped)
    overlapped_df[, c(1, 2, 3, 4, 5, 6)] <- lapply(overlapped_df[, c(1, 2, 3, 4, 5, 6)], as.character)
    # form full list of multiple overlapping
    overlap_temp <- list()
    for (x in seq(1, nrow(overlapped_df), 1)) {
      df_temp <- overlapped_df[x, c("seqnames", "start", "end", "names")]
      for (y in overlapped_df[x, c("motif_id")][[1]]){
        df_temp$motif_id <- y
        overlap_temp[[paste0(x, "_", y)]] <- df_temp
      }
    }
    overlapped_df_full <- as.data.frame(matrix(nrow = length(overlap_temp), ncol = 5))
    for (x in seq(1, nrow(overlapped_df_full), 1)) {
      overlapped_df_full[x, ] <- overlap_temp[[x]]
    }
    colnames(overlapped_df_full) <- c(
      "peak_chr","peak_start","peak_end","peak_id","motif_id"
    )
    overlapped_df_full_all <- merge(
      x = overlapped_df_full,
      y = motif_seq[, c("motif_chr", "motif_start", "motif_end", "motif_id")],
      by = "motif_id"
    )
    overlapped_df_full_all[, c("peak_start", "peak_end", "motif_start", "motif_end")] <- lapply(overlapped_df_full_all[, c("peak_start","peak_end","motif_start","motif_end")], as.numeric)
    overlapped_df_full_all$peak_center <- (overlapped_df_full_all$peak_start + overlapped_df_full_all$peak_end)/2
    overlapped_df_full_all[, c("peak_center")] <- vapply(overlapped_df_full_all[,c("peak_center")], as.integer, integer(1))
    overlapped_df_full_all$left_point <- overlapped_df_full_all$motif_start-overlapped_df_full_all$peak_center
    overlapped_df_full_all$right_point <- overlapped_df_full_all$motif_end-overlapped_df_full_all$peak_center
    #calculate motif occurrence
    motif_occurrence <- c(rep(0,(as.integer(plot_at_each_side)*2+3)))
    for (x in seq(1, nrow(overlapped_df_full_all), 1)){
      left_point <- overlapped_df_full_all[x, c("left_point")]
      right_point <- overlapped_df_full_all[x, c("right_point")]
      left_point <- left_point + (as.integer(plot_at_each_side)+2)
      left_point <- max(1, left_point)
      left_point <- min((as.integer(plot_at_each_side)*2+3), left_point)
      right_point <- right_point + (as.integer(plot_at_each_side)+2)
      right_point <- max(1, right_point)
      right_point <- min((as.integer(plot_at_each_side)*2+3), right_point)
      motif_occurrence[left_point:right_point] <- motif_occurrence[left_point:right_point]+1
    }
    motif_occurrence <- 100*motif_occurrence/nrow(overlapped_df_full_all)
    motif_occurrence <- motif_occurrence[2:(as.integer(plot_at_each_side)*2+2)]
    motifDistrb_list_add <- list("target_id"=id,"input_peak_number"=nrow(peak_target),"peak_with_motif"=nrow(overlapped_df), "motif_occurrence"=motif_occurrence)
    motifDistrb_list[[peak_id_i]] <- motifDistrb_list_add
  }
  return(motifDistrb_list)
}
