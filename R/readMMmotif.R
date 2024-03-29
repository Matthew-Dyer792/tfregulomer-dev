readMMmotif <- function(motif_file_path, motif_format, id, num_peak,
                        nsites, motif_file_path_MEME)
{
  con <- tryCatch({
    file(motif_file_path, "r")
  },
  warning = function(w) {
    message(paste0("Warning: Please check your id, motif format or internet access!"))
    stop()
  },
  error = function(err) {
    message(paste0("Warning: Please check your id, motif format or internet access!"))
  })

  # start to process motif file if any
  if (!is.null(con)) {
    # count motif file line
    line_count <- 0
    # store motif matrix values
    motif_vector <- c()
    if (motif_format == "MEME") {
      lines <- readLines(con)
      if (length(lines) != 0) {
        # version line
        line_vector <- unlist(strsplit(as.character(lines[[1]]), split = " "))
        line_vector <- line_vector[line_vector != ""]
        version <- as.numeric(substr(line_vector[[3]], start = 1, stop = 1))

        # alphabet line
        line_vector <- unlist(strsplit(as.character(lines[[3]]), split = "="))
        line_vector <- line_vector[line_vector != ""]
        alphabet <- as.character(trimws(
          tail(line_vector, n = 1), which = "both"
        ))

        # strands line
        line_vector <- unlist(strsplit(as.character(lines[[5]]), split = ":"))
        line_vector <- line_vector[line_vector != ""]
        strand <- as.character(trimws(tail(line_vector, n = 1), which = "both"))

        # background line
        line_vector <- unlist(strsplit(as.character(lines[[8]]), split = " "))
        line_vector <- line_vector[line_vector != ""]
        background <- c(
          "A" = as.numeric(line_vector[2]), "C" = as.numeric(line_vector[4]),
          "G" = as.numeric(line_vector[6]), "T" = as.numeric(line_vector[8])
        )

        # id line
        line_vector <- unlist(strsplit(as.character(lines[[10]]), split = " "))
        line_vector <- line_vector[line_vector != ""]
        if (nchar(line_vector[2]) > nchar(line_vector[3])) {
          MMmotif_id <- id
          alternate_name <- as.character(line_vector[3])
        } else {
          MMmotif_id <- id
          alternate_name <- as.character(line_vector[2])
        }

        # width, evalue
        line_vector <- unlist(strsplit(as.character(lines[[12]]), split = " "))
        line_vector <- line_vector[line_vector != ""]
        width <- as.integer(line_vector[6])
        evalue <- as.numeric(line_vector[10])

        # motif matrix
        for (line in lines[-1:-12]) {
          line_vector <- unlist(strsplit(as.character(line), split = " "))
          line_vector <- line_vector[line_vector != ""]
          for (line_vector_item in line_vector) {
            motif_vector <- c(motif_vector, as.numeric(trimws(
              line_vector_item, which = "both"
            )))
          }
        }
      }

      motif_matrix <- matrix(motif_vector, ncol = 4, byrow = TRUE)
      colnames(motif_matrix) <- c("A", "C", "G", "T")
      close(con)
      MMmotif <- new("MMmotif")
      MMmotif <- updateMMmotif(MMmotif, motif_format = motif_format,
                               version = version, alphabet = alphabet,
                              strand = strand, background = background,
                              id = MMmotif_id, alternate_name = alternate_name,
                              width = width, nsites = nsites, nPeaks = num_peak,
                              evalue = evalue, motif_matrix = motif_matrix )
    } else {
      while (TRUE) {
        line <- readLines(con, n=1)
        line_count <- line_count + 1
        if (line_count > 6 && grepl(pattern = "XX", line)) {
          break
        }
        # id line
        if (line_count == 5) {
          line_vector <- unlist(strsplit(as.character(line), split=" "))
          line_vector <- line_vector[line_vector != ""]
          if (nchar(line_vector[2]) > nchar(line_vector[3])) {
            MMmotif_id <- id
            alternate_name <- as.character(line_vector[3])
          } else {
            MMmotif_id <- id
            alternate_name <- as.character(line_vector[2])
          }
        }
        # motif matrix
        if (line_count > 6 & !grepl(pattern = "XX", line)) {
          line_vector <- unlist(strsplit(as.character(line), split = "\t"))
          line_vector <- line_vector[line_vector != ""]
          for (line_vector_item in line_vector) {
            motif_vector <- c(motif_vector, as.numeric(trimws(line_vector_item, which = "both")))
          }
        }
      }
      motif_matrix <- matrix(motif_vector, ncol = 5, byrow = TRUE)
      motif_matrix <- motif_matrix[,-1]
      colnames(motif_matrix) <- c("A", "C", "G", "T")
      width <- as.integer(nrow(motif_matrix))
      close(con)

      # get background from MEME
      con_MEME <- file(motif_file_path_MEME, "r")
      line_count_MEME <- 0
      background <- c("A" = 0.25, "C" = 0.25, "G" = 0.25, T = 0.25)
      while (TRUE) {
        line <- readLines(con, n = 1)
        line_count_MEME <- line_count_MEME + 1
        if (length(line) == 0) {
          break
        }
        # background line
        if (line_count_MEME == 8) {
          line_vector <- unlist(strsplit(as.character(line), split = " "))
          line_vector <- line_vector[line_vector != ""]
          background <- c(
            "A" = as.numeric(line_vector[2]), "C" = as.numeric(line_vector[4]),
            "G" = as.numeric(line_vector[6]), "T" = as.numeric(line_vector[8])
          )
        }
        # width, evalue
        if (line_count_MEME == 12) {
          line_vector <- unlist(strsplit(as.character(line), split = " "))
          line_vector <- line_vector[line_vector != ""]
          evalue <- as.numeric(line_vector[10])
          break
        }
      }
      close(con_MEME)

      MMmotif <- new("MMmotif")
      MMmotif <- updateMMmotif(MMmotif, motif_format = motif_format,
                               version = 0, alphabet = "ACGT",
                              strand = "+ -", background = background,
                              id = MMmotif_id, alternate_name = alternate_name,
                              width = width, nsites = nsites, nPeaks = num_peak,
                              evalue = evalue, motif_matrix = motif_matrix)
    }

    return(MMmotif)
  } else {
    message("The motif file is NULL!")
    return(NULL)
  }
}
