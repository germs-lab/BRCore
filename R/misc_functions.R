# Miscellaneous functions for Inter-BRC Core Microbiome Analysis
# By Bolívar Aponte Rolón

# extract_matrix(): extract/subset a ASV/OTU matrix based on a vector of strings and a rowSums() condition.

extract_matrix <- function(physeq, .vec, keep_rows_sums = 0) {
  # Error handling
  if (!inherits(physeq, c("phyloseq", "matrix"))) {
    cli::cli_abort(
      "{.arg physeq} must be a 'phyloseq' or 'matrix' object.\nYou've supplied a {class(physeq)[1]} class."
    )
  }

  if (inherits(physeq, "phyloseq")) {
    cli::cli_alert_info("Detected a 'phyloseq' object. Input object is valid!")
    physeq <- otu_table(physeq)
  } else {
    cli::cli_alert_info("Detected a 'matrix' object. Proceeding with the input matrix.")
  }

  if (!inherits(physeq, "matrix")) {
    cli::cli_abort(
      "The extracted OTU/ASV table is not a matrix. Something went wrong."
    )
  }

  cli::cli_alert_success("Extracting the ASV/OTU matrix...")

  # Main function
  physeq %>%
    t() %>% # We need samples as rows and ASV as columns
    as.data.frame() %>%
    select(contains(.vec)) %>%
    .[rowSums(.) > keep_rows_sums, ] %>% # Keep only samples with a non-zero sum. Not all samples have the "core".
    as.matrix()

  # cli::cli_alert_success("Yay! You got a matrix.")
}



# Borrowed subset.fasta from https://github.com/GuillemSalazar/FastaUtils/blob/master/R/FastaUtils_functions.R

#' Select a subset of sequences from a fasta file
#'
#' This function loads a fasta file, selects a subset on sequences based on the sequence's names and writes a new fsata file.
#' @param file Fasta file
#' @param subset Vector with (exact) names of the sequences to be retrieved.
#' @param out Path to output file. If absent, the '.mod' is added to the input file's name (path is thus conserved).
#' @keywords FastaUtils
#' @return Writes a fasta file with the selected sequences.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' subset_fasta(
#' file = "http://greengenes.lbl.gov/Data/JD_Tutorial/UnAligSeq24606.txt", 
#' subset = c("24.6jsd1.Tut", "24.6jsd2.Tut ", "24.6jsd3.Tut "), 
#' out = "out.fasta"
#' )
subset_fasta <- function(file = NULL,
                         subset = NULL,
                         out = paste(file, ".subset", sep = "")) {
  library(Biostrings)
  sequences <- readDNAStringSet(file)
  if (all(as.character(subset) %in% names(sequences)) == FALSE)
    stop("There are names in 'subset' not present in the fasta file")
  pos <- match(as.character(subset), names(sequences))
  writeXStringSet(sequences[pos], filepath = out)
}


# Parallelized function for choosing dimensions for NMDS
nmds_screen_parallel <- function(x, ncores = parallel::detectCores() - 1) {
  # Function to calculate stress for a given number of dimensions
  calculate_stress <- function(k) {
    replicate(10, metaMDS(x, autotransform = FALSE, k = k, maxit = 100, trymax = 10)$stress)
  }

  # Use mclapply to parallelize the stress calculation
  stress_values <- parallel::mclapply(1:10, calculate_stress, mc.cores = ncores)

  # Plot the results
  plot(rep(1, 10), stress_values[[1]],
    xlim = c(1, 10),
    ylim = c(0, 0.30),
    xlab = "# of Dimensions",
    ylab = "Stress",
    main = "NMDS stress plot"
  )

  for (i in 1:9) {
    points(rep(i + 1, 10), stress_values[[i + 1]])
  }
}

# Standard NMDS plots
gg_nmds <- function(.data, .color, .shape = NULL, .drop_na) {
  .color <- enquo(.color)
  .drop_na <- enquo(.drop_na)

  .data %>%
    drop_na(!!.drop_na) %>%
    ggplot(., aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = !!.color, stroke = 1),
      alpha = 0.5,
      na.rm = TRUE
    ) +
    geom_hline(
      yintercept = 0,
      colour = "grey70",
      linewidth = 0.65
    ) +
    geom_vline(
      xintercept = 0,
      colour = "grey70",
      linewidth = 0.65
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text()
    ) +
    stat_ellipse(
      aes(color = !!.color),
      geom = "path",
      linewidth = 1.3,
      position = "identity",
      type = "t",
      linetype = 1,
      level = 0.95,
      segments = 51,
      na.rm = FALSE,
      show.legend = NA,
      inherit.aes = TRUE
    )
}
