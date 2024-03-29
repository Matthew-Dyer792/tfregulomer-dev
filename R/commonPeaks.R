#' commonPeaks
#'
#' This function allows you to obtain a list of common peak subsets along with the DNA methylation profiles.
#' @param target_peak_id Character of vector, each of which is a TFregulomeR ID. Each of target peak will be compared with all "compared peaks" to get its common subset.
#' @param motif_only_for_target_peak Either TRUE of FALSE (default). If TRUE, only peaks with motif will be loaded for each TFregulomeR ID in target_peak_id.
#' @param user_target_peak_list A list of data.frames, each of which contains user's own bed-format target peak regions.
#' @param user_target_peak_id Character of vector, each of which is a unique ID corresponding to each peak set in the list user_target_peak_list. If the IDs are not provided or not unique, the function will automatically generate the IDs of its own. If any of the peak sets is derived from TFregulomeR, its TFregulomeR ID should be used here correspondingly.
#' @param compared_peak_id Character of vector, each of which is a TFregulomeR ID.
#' @param motif_only_for_compared_peak Either TRUE of FALSE (default). If TRUE, only peaks with motif will be loaded for each TFregulomeR ID in compared_peak_id.
#' @param user_compared_peak_list A list of data.frames, each of which contains user's own bed-format compared peak regions.
#' @param user_compared_peak_id Character of vector, each of which is a unique ID corresponding to each peak set in the list user_compared_peak_list. If the IDs are not provided or not unique, the function will automatically generate the IDs of its own. If any of the peak sets is derived from TFregulomeR, its TFregulomeR ID should be used here correspondingly.
#' @param methylation_profile_in_narrow_region Either TRUE (default) of FALSE. If TRUE, methylation states in 200bp window surrounding peak summits for each common peak from target_peak_id and user_target_peak_list with TFregulomeR ID.
#' @param motif_type Motif PFM format, either in MEME by default or TRANSFAC.
#' @param server server localtion to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "https://bioinfo-csi.nus.edu.sg/methmotif/" or "https://methmotif.org", please use a new url.
#' @param local_db_path The complete path to the SQLite implementation of TFregulomeR database available at "https://methmotif.org/API_ZIPPED.zip"
#' @return  matrix of CommonPeaksMM class objects
#' @keywords commonPeaks
#' @export
#' @examples
#' target_id <- c("MM1_HSA_K562_CEBPB")
#' compared_id <- c("MM1_HSA_HepG2_CEBPB")
#' commonPeaks_output <- commonPeaks(target_peak_id=target_id,
#'                                   motif_only_for_target_peak=TRUE,
#'                                   compared_peak_id=compared_id,
#'                                   motif_only_for_compared_peak=TRUE,
#'                                   methylation_profile_in_narrow_region=TRUE)

commonPeaks <- function(target_peak_id,
                        motif_only_for_target_peak = FALSE,
                        user_target_peak_list,
                        user_target_peak_id,
                        compared_peak_id,
                        motif_only_for_compared_peak = FALSE,
                        user_compared_peak_list,
                        user_compared_peak_id,
                        methylation_profile_in_narrow_region = TRUE,
                        motif_type = "MEME",
                        server = "ca",
                        TFregulome_url,
                        local_db_path = NULL)
{
  # check the input arguments
  if (missing(target_peak_id) && missing(user_target_peak_list)) {
    stop("No target peak input. Please input TFregulomeR peaks using TFregulomeR ID(s) by 'target_peak_id = ' OR your own peak list using a list of data.frame(s) containing bed-format regions by 'user_target_peak_list = '")
  }
  if (missing(compared_peak_id) && missing(user_compared_peak_list)) {
    stop("No compared peak input. Please input TFregulomeR peaks using TFregulomeR ID(s) by 'compared_peak_id = ' OR your own peak list using a list of data.frame(s) containing bed-format regions by 'user_compared_peak_list = '")
  }
  if ((!missing(user_target_peak_list) && !is.list(user_target_peak_list)) ||
      (!missing(user_compared_peak_list) && !is.list(user_compared_peak_list))) {
    stop("The class of input 'user_target_peak_list' and 'user_compared_peak_list' should be 'list', a list of bed-like data.frame storing peak regions!")
  }
  if (!is.logical(motif_only_for_target_peak) || !is.logical(motif_only_for_compared_peak)) {
   stop("motif_only_for_target_peak and motif_only_for_compared_peak should be either TRUE or FALSE (default)")
  }
  if (!is.logical(methylation_profile_in_narrow_region)) {
    stop("methylation_profile_in_narrow_region should be either TRUE (default) or FALSE")
  }
  if (motif_type != "MEME" && motif_type != "TRANSFAC") {
    stop("motif_type should be either 'MEME' (default) or 'TRANSFAC'!")
  }

  # build api_object
  api_object <- .construct_api(server, TFregulome_url, local_db_path)

  message("TFregulomeR::commonPeaks() starting ... ...")
  if (methylation_profile_in_narrow_region) {
    message("You chose to profile the methylation levels in 200bp window around peak summits, if there is any peak set loaded from TFregulomeR")
  } else {
    message("You chose NOT to profile the methylation levels in 200bp window around peak summits")
  }
  # loading target peak list
  message("Loading target peak list ... ...")
  target_peak_list_all <- list()
  # loading from TFregulomeR server
  TFregulome_target_peak_id <- c()
  is_taregt_TFregulome <- c()
  target_list_count <- 0
  if (!missing(target_peak_id) && length(target_peak_id) > 0) {
    message(paste0("... You have ", length(target_peak_id)," TFBS(s) requested to be loaded from TFregulomeR server"))
    if (motif_only_for_target_peak == TRUE) {
      message("... You chose to load TF peaks with motif only. Using 'motif_only_for_target_peak' tunes your options")
    } else {
      message("... You chose to load TF peaks regardless of presence of motif. Using 'motif_only_for_target_peak' tunes your options")
    }
    message("... loading TFBS(s) from TFregulomeR now")
    for (i in target_peak_id) {
      peak_i <- suppressMessages(.loadPeaks(
        id = i,
        includeMotifOnly = motif_only_for_target_peak,
        api_object = api_object
      ))
      if (is.null(peak_i)) {
        message(paste0("... ... NO peak file for your id '", i,"'."))
      } else {
        target_list_count <- target_list_count + 1
        target_peak_list_all[[target_list_count]] <- peak_i
        TFregulome_target_peak_id <- c(TFregulome_target_peak_id, i)
        is_taregt_TFregulome <- c(is_taregt_TFregulome, TRUE)
        message(paste0("... ... peak file loaded successfully for id '", i,"'"))
      }
    }
    message("... Done loading TFBS(s) from TFregulomeR")
  }
  # users' peaks
  if (!missing(user_target_peak_list) && length(user_target_peak_list) > 0) {
    message(paste0("... You have ",length(user_target_peak_list)," customised peak set(s)"))
    if (missing(user_target_peak_id) || length(user_target_peak_id) != length(user_target_peak_list) ||
        length(unique(user_target_peak_id)) != length(user_target_peak_list)) {
      message("... ... You didn't provide the ID for each customised peak set or your ID number does not uniquely equal to the input user peak number. Instead we will use 'user_target_peak1', 'user_target_peak2'..." )
      user_target_peak_id <- paste0("user_target_peak", seq(1,length(user_target_peak_list), 1))
    }
    # new user target peak id in case that any input peak set is empty
    user_target_peak_id_new <- c()
    for (i in seq(1, length(user_target_peak_list), 1)) {
      peak_i <- user_target_peak_list[[i]]
      if (nrow(peak_i) == 0) {
        message(paste0("... ... Your input peak set '",user_target_peak_id[i],"' is empty, so SKIP!"))
      } else {
        user_target_peak_id_new <- c(user_target_peak_id_new, user_target_peak_id[i])
        colname_new <- colnames(peak_i)
        colname_new[1] <- "chr"
        colname_new[2] <- "start"
        colname_new[3] <- "end"
        # check if the input peak set has fourth column for id
        no_id <- TRUE
        if (length(colname_new) >= 4) {
          if (length(unique(peak_i[, 4])) == nrow(peak_i)) {
            colname_new[4] <- "id"
            no_id <- FALSE
          }
        }
        colnames(peak_i) <- colname_new
        if (no_id) {
          peak_i$id <- paste0(
            user_target_peak_id[i],
            "_",
            as.vector(rownames(peak_i))
          )
        }
        target_list_count <- target_list_count + 1
        target_peak_list_all[[target_list_count]] <- peak_i
        # test if user input id i match any TFregulomeR ID
        motif_matrix_i <- suppressMessages(.searchMotif(
          id = user_target_peak_id[i],
          motif_format = motif_type,
          api_object = api_object
        ))
        if (is.null(motif_matrix_i)) {
          is_taregt_TFregulome <- c(is_taregt_TFregulome, FALSE)
        } else {
          is_taregt_TFregulome <- c(is_taregt_TFregulome, TRUE)
        }
      }
    }
  } else {
    user_target_peak_id_new <- c()
  }
  # combine TFregulomeR ID and user ID
  target_peak_id_all <- c(TFregulome_target_peak_id, user_target_peak_id_new)

  # loading compared peak list
  message("Loading compared peak list ... ...")
  compared_peak_list_all <- list()
  is_compared_TFregulome <- c()
  # loading from TFregulomeR server
  compared_list_count <- 0
  TFregulome_compared_peak_id <- c()
  if (!missing(compared_peak_id) && length(compared_peak_id) > 0) {
    message(paste0("... You have ", length(compared_peak_id)," TFBS(s) requested to be loaded from TFregulomeR server"))
    if (motif_only_for_compared_peak == TRUE) {
      message("... You chose to load TF peaks with motif only. Using 'motif_only_for_compared_peak' tunes your options")
    } else {
      message("... You chose to load TF peaks regardless of presence of motif. Using 'motif_only_for_compared_peak' tunes your options")
    }
    message("... loading TFBS(s) from TFregulomeR now")
    for (i in compared_peak_id) {
      peak_i <- suppressMessages(.loadPeaks(
        id = i,
        includeMotifOnly = motif_only_for_compared_peak,
        api_object = api_object
      ))
      if (is.null(peak_i)) {
        message(paste0("... ... NO peak file for your id '", i,"'."))
      } else {
        compared_list_count <- compared_list_count + 1
        compared_peak_list_all[[compared_list_count]] <- peak_i
        is_compared_TFregulome <- c(is_compared_TFregulome, TRUE)
        TFregulome_compared_peak_id <- c(TFregulome_compared_peak_id, i)
        message(paste0("... ... peak file loaded successfully for id '", i,"'"))
      }
    }
    message("... Done loading TFBS(s) from TFregulomeR")
  }
  # users' peaks
  if (!missing(user_compared_peak_list) && length(user_compared_peak_list) > 0) {
    message(paste0("... You have ",length(user_compared_peak_list)," customised peak set(s)"))
    if (missing(user_compared_peak_id) || length(user_compared_peak_id)!=length(user_compared_peak_list) ||
        length(unique(user_compared_peak_id))!=length(user_compared_peak_list))
    {
      message("... ... You didn't provide the ID for each customised peak set or your ID number does not uniquely equal to the input user peak number. Instead we will use 'user_compared_peak1', 'user_compared_peak2'..." )
      user_compared_peak_id <- paste0("user_compared_peak", seq(1,length(user_compared_peak_list), 1))
    }
    # new user compared peak id in case that any input peak set is empty
    user_compared_peak_id_new <- c()
    for (i in seq(1, length(user_compared_peak_list), 1)) {
      peak_i <- user_compared_peak_list[[i]]
      if (nrow(peak_i) == 0) {
        message(paste0("... ... Your input peak set '",user_compared_peak_id[i],"' is empty, so SKIP!"))
      } else {
        user_compared_peak_id_new <- c(user_compared_peak_id_new, user_compared_peak_id[i])
        peak_i_sub <- peak_i[,c(1,2,3)]
        colnames(peak_i_sub) <- c("chr","start","end")
        peak_i_sub$id <- paste0("compared_peak_", as.vector(rownames(peak_i_sub)))
        peak_i_sub <- peak_i_sub[,c("chr","start","end","id")]
        compared_list_count <- compared_list_count + 1
        compared_peak_list_all[[compared_list_count]] <- peak_i_sub
        # test if user input id i match any TFregulomeR ID
        motif_matrix_i <- suppressMessages(.searchMotif(
          id = user_compared_peak_id[i],
          motif_format = motif_type,
          api_object = api_object
        ))
        if (is.null(motif_matrix_i)) {
          is_compared_TFregulome <- c(is_compared_TFregulome, FALSE)
        } else {
          is_compared_TFregulome <- c(is_compared_TFregulome, TRUE)
        }
      }
    }
  }

  # check if target peak list is empty
  if (length(target_peak_list_all) == 0) {
    message("No input in the target peak list. Function ends here!")
    return(NULL)
  }
  # check if compared peak list is empty
  if (length(compared_peak_list_all) == 0) {
    message("No input in the compared peak list. Function ends here!")
    return(NULL)
  }
  # start analysing
  common_peak_matrix <- list()
  for (i in seq(1, length(target_peak_list_all), 1)) {
    target_id_i <- target_peak_id_all[i]
    target_peak_i <- target_peak_list_all[[i]]
    number_of_orignal_target <- nrow(target_peak_i)
    message(paste0("Start analysing: ", target_id_i, "... ..."))

    ## if it is from TFregulomeR
    if (is_taregt_TFregulome[i]) {
      isTFregulome_target <- TRUE
      # make the request
      request_content_df <- apiRequest(api_object, id = target_id_i)
      source_i <- request_content_df[, "source"]
      if (source_i == "MethMotif") {
        isMethMotifID_target <- TRUE
        motif_seq_path_target <- request_content_df[1, c("TFBS")]
        meth_file_path_target <- request_content_df[1, c("DNA_methylation_profile")]
        meth_file_200bp_path_target <- request_content_df[1, c("DNA_methylation_profile_200bp")]
        WGBS_replicate_target <- request_content_df[1, c("WGBS_num")]
      } else {
        isMethMotifID_target <- FALSE
        motif_seq_path_target <- request_content_df[1, c("TFBS")]
      }
    } else {
      isTFregulome_target <- FALSE
    }
    # comparing with compared peak list
    for (j in seq(1, length(compared_peak_list_all), 1)) {
      compared_peak_j <- compared_peak_list_all[[j]]
      if (isTFregulome_target) {
        bed_target_i <- GRanges(target_peak_i$chr,
                                IRanges(target_peak_i$start-99,
                                        target_peak_i$end+100),
                                id=target_peak_i$id)
      } else {
        bed_target_i <- GRanges(target_peak_i$chr,
                                IRanges(target_peak_i$start,
                                        target_peak_i$end),
                                id=target_peak_i$id)
      }
      if (is_compared_TFregulome[j])
      {
        bed_compared_j <- GRanges(compared_peak_j$chr,
                                  IRanges(compared_peak_j$start-99,
                                          compared_peak_j$end+100),
                                  id=compared_peak_j$id)
      }
      else
      {
        bed_compared_j <- GRanges(compared_peak_j$chr,
                                  IRanges(compared_peak_j$start,
                                          compared_peak_j$end),
                                  id=compared_peak_j$id)
      }
      # get target peak intersecting with compared peak
      # subsetOverlaps may mis-think the two sets coming from different references, so suppressWarnings here
      suppressWarnings(bedTarget_with_bedcompared <- subsetByOverlaps(bed_target_i, bed_compared_j))
      peakTarget_with_peakcompared <- unique(as.data.frame(bedTarget_with_bedcompared))
      target_peak_with_peakcompared <- target_peak_i[which(target_peak_i$id %in% peakTarget_with_peakcompared$id), ]
      target_peak_i <- target_peak_with_peakcompared
    }

    MethMotif_target <- new("MethMotif")
    #initiate methylation profile distribution
    meth_score_distri_target <- matrix()

    if (isTFregulome_target)
    {
      motif_seq_target <- read.delim(motif_seq_path_target, sep = "\t", header = FALSE)
      if (nrow(target_peak_i)>0)
      {
        #compute motif matrix
        colnames(motif_seq_target) <- c("chr","start","end","strand","weight", "pvalue","qvalue","sequence")
        motif_len_target <- nchar(as.character(motif_seq_target[1,"sequence"]))
        motif_seq_target$id <- paste0(target_id_i,"_motif_sequence_", as.vector(rownames(motif_seq_target)))
        motif_seq_target_grange <- GRanges(motif_seq_target$chr,
                                           IRanges(motif_seq_target$start+1,
                                                   motif_seq_target$end),
                                           id=motif_seq_target$id)
        bed_target_done_common <- GRanges(target_peak_i$chr,
                                          IRanges(target_peak_i$start-99,
                                                  target_peak_i$end+100),
                                          id=target_peak_i$id)
        suppressWarnings(motif_of_bed_target_done_common <- subsetByOverlaps(motif_seq_target_grange, bed_target_done_common))
        motif_of_peakTarget_done_common <- unique(as.data.frame(motif_of_bed_target_done_common))

        if (nrow(motif_of_peakTarget_done_common)>0)
        {
          # nPeaks
          suppressWarnings(peaks_with_motif_gr <- subsetByOverlaps(bed_target_done_common,
                                                                   motif_seq_target_grange))
          peaks_with_motif_df <- unique(as.data.frame(peaks_with_motif_gr))
          number_of_target_peak_i_with_motif <- nrow(peaks_with_motif_df)


          motif_of_peakTarget_done_common_allInfo <- motif_seq_target[which(motif_seq_target$id %in% motif_of_peakTarget_done_common$id),]
          motif_matrix_of_peakTarget_done_common <- formMatrixFromSeq(input_sequence = as.vector(motif_of_peakTarget_done_common_allInfo$sequence),
                                                                      motif_format = motif_type)
          # compute beta score matrix
          if (isMethMotifID_target)
          {
            # calculate beta score matrix
            # methylation file can be empty
            meth_level_target <- tryCatch(read.delim(meth_file_path_target, sep = "\t", header = FALSE),
                                          error=function(e) data.frame())
            # methylation file can be empty
            if (nrow(meth_level_target)==0)
            {
              beta_score_matrix_of_peakTarget_done_common <- formBetaScoreFromSeq(input_meth = data.frame(),
                                                                                  WGBS_replicate = WGBS_replicate_target,
                                                                                  motif_len = motif_len_target)
            }
            else
            {
              colnames(meth_level_target) <- c("chr","start","end","meth_score","C_num","T_num","seq_chr","seq_start",
                                               "seq_end","strand","weight","pvalue","qvalue","sequence")
              meth_level_target$id <- paste0(target_id_i,"_motif_with_CG_", as.vector(rownames(meth_level_target)))
              meth_level_target_grange <- GRanges(meth_level_target$seq_chr,
                                                  IRanges(meth_level_target$seq_start,
                                                          meth_level_target$seq_end),
                                                  id=meth_level_target$id)
              suppressWarnings(meth_level_peakTarget_done_common <- unique(as.data.frame(subsetByOverlaps(meth_level_target_grange,
                                                                         motif_of_bed_target_done_common))))
              meth_level_peakTarget_done_common_allInfo <- meth_level_target[which(meth_level_target$id %in% meth_level_peakTarget_done_common$id),]
              beta_score_matrix_of_peakTarget_done_common <- formBetaScoreFromSeq(input_meth = meth_level_peakTarget_done_common_allInfo,
                                                                                  WGBS_replicate = WGBS_replicate_target,
                                                                                  motif_len = motif_len_target)
            }
          }
          else
          {
            beta_score_matrix_of_peakTarget_done_common <- as.matrix(NA)
          }

          if (motif_type == "TRANSFAC")
          {
            version_target <- 0
          }
          else
          {
            version_target <- 4
          }
          MethMotif_target@MMmotif <- updateMMmotif(MethMotif_target@MMmotif,
                                                    motif_format = motif_type,
                                                    version = version_target,
                                                    background = c("A"=0.25,"C"=0.25,"G"=0.25,"T"=0.25),
                                                    id = paste0(target_id_i,"_common_peaks"),
                                                    alternate_name = target_id_i,
                                                    width = motif_len_target,
                                                    nsites = nrow(motif_of_peakTarget_done_common),
                                                    nPeaks = number_of_target_peak_i_with_motif,
                                                    motif_matrix=motif_matrix_of_peakTarget_done_common)
          MethMotif_target@MMBetaScore <- beta_score_matrix_of_peakTarget_done_common
        }


        # profile methylation level in narrow regions
        if (methylation_profile_in_narrow_region)
        {
          ### if in 200bp around peaks
          if (isMethMotifID_target)
          {
            meth_level_200bp_target <- tryCatch(read.delim(meth_file_200bp_path_target, sep = "\t", header = FALSE),
                                                error=function(e) data.frame())
            if (nrow(meth_level_200bp_target) == 0)
            {
              meth_score_distri_target <- formBetaScoreDistri(input_meth = data.frame())
            }
            else
            {
              colnames(meth_level_200bp_target) <- c("chr","start","end",
                                                     "meth_score","C_num","T_num")
              meth_level_200bp_target$id <- paste0(target_id_i,"_200bp_CG_", as.vector(rownames(meth_level_200bp_target)))
              meth_level_200bp_target_grange <- GRanges(meth_level_200bp_target$chr,
                                                        IRanges(meth_level_200bp_target$start,
                                                                meth_level_200bp_target$end),
                                                        id=meth_level_200bp_target$id)
              bed_target_i <- GRanges(target_peak_i$chr,
                                      IRanges(target_peak_i$start-99,
                                              target_peak_i$end+100),
                                      id=target_peak_i$id)
              suppressWarnings(meth_level_in_common_peaks_200bp <- unique(as.data.frame(subsetByOverlaps(meth_level_200bp_target_grange,
                                                                                                         bed_target_i))))
              meth_level_in_common_peaks_200bp_allInfo <- unique(meth_level_200bp_target[which(meth_level_200bp_target$id
                                                                                               %in% meth_level_in_common_peaks_200bp$id),])
              meth_score_distri_target <- formBetaScoreDistri(input_meth = as.data.frame(meth_level_in_common_peaks_200bp_allInfo$meth_score))
            }
          }
        }
      }
    }
    new_CommonPeaksMM <- new("CommonPeaksMM")
    new_CommonPeaksMM <- updateCommonPeaksMM(theObject = new_CommonPeaksMM,
                                             id = paste0(target_id_i, "_common_peaks"),
                                             common_percentage = 100*nrow(target_peak_i)/number_of_orignal_target,
                                             common_peak = target_peak_i,
                                             isTFregulomeID = isTFregulome_target,
                                             MethMotif = MethMotif_target,
                                             methylation_profile = meth_score_distri_target)
    common_peak_matrix[[paste0(target_id_i, "_common_peaks")]] <- new_CommonPeaksMM
  }
  message("Done analysing.")
  dim(common_peak_matrix) <- c(length(target_peak_id_all), 1)
  rownames(common_peak_matrix) <- c(target_peak_id_all)
  colnames(common_peak_matrix) <- c("common")
  return(common_peak_matrix)
}


formBetaScoreDistri <- function(input_meth)
{
  beta_score_distri <- c(rep(0,10))
  names(beta_score_distri) <- c("0-10%","10-20%","20-30%","30-40%","40-50%",
                                "50-60%","60-70%","70-80%","80-90%","90-100%")
  if (nrow(input_meth)>0)
  {
    colnames(input_meth) <- c("beta_score")
    betascore_hist <- hist(input_meth$beta_score, breaks = c(seq(0,100,10)), plot = FALSE)
    beta_score_distri <- c(betascore_hist$counts)
    names(beta_score_distri) <- c("0-10%","10-20%","20-30%","30-40%","40-50%",
                                  "50-60%","60-70%","70-80%","80-90%","90-100%")
  }
  beta_score_distri_matrix <- as.matrix(beta_score_distri)
  colnames(beta_score_distri_matrix) <- "CpG_num"
  return(beta_score_distri_matrix)
}




