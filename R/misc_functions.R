#' Extract ASV/OTU Matrix from phyloseq or Matrix Object
#'
#' This function extracts an ASV/OTU matrix from a `phyloseq` or `matrix` object. 
#' It ensures the matrix is in the correct format (samples as rows and ASVs as columns) 
#' and filters samples based on a specified vector and row sums.
#'
#' @param physeq A `phyloseq` or `matrix` object containing ASV/OTU data. 
#'   If a `phyloseq` object is provided, the OTU table is extracted. 
#'   If a `matrix` is provided, it is used directly.
#' @param .vec A character string or vector specifying the columns to select from the matrix. 
#'   This is passed to `dplyr::select(contains(.vec))` to filter the matrix.
#' @param keep_rows_sums A numeric value specifying the minimum row sum required to keep a sample. 
#'   Rows (samples) with sums less than or equal to this value are removed. Default is `0`.
#' @param taxa_are_rows A logical value specifying whether the taxa are rows in the matrix. Default is `TRUE`.
#' @return A matrix with samples as rows and ASVs/OTUs as columns. Only samples matching 
#'   the `.vec` criteria and with row sums greater than `keep_rows_sums` are retained.
#'
#' @import phyloseq
#' @importFrom dplyr select contains
#' @importFrom cli cli_abort cli_alert_info cli_alert_success
#' @export
#'
#' @examples
#' # Load phyloseq object
#' data(GlobalPatterns, package = "phyloseq")
#'
#' # Extract matrix for samples containing "Feces" in their names
#' fec_matrix <- extract_matrix(GlobalPatterns, .vec = "Feces", keep_rows_sums = 10)
#'
#' # View the extracted matrix
#' head(fec_matrix)

extract_matrix <- function(physeq, .vec, keep_rows_sums = 0, taxa_are_rows = TRUE) {
  # Error handling
  if (!inherits(physeq, c("phyloseq", "matrix"))) {
    cli::cli_abort(
      "{.arg physeq} must be a 'phyloseq' or 'matrix' object.\nYou've supplied a {class(physeq)[1]} class."
    )
  }

  if (inherits(physeq, "phyloseq")) {
    cli::cli_alert_info("Detected a 'phyloseq' object. Input object is valid!")
    physeq <- otu_table(physeq, taxa_are_rows = taxa_are_rows)
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
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>

subset_fasta <- function(file = NULL,
                         subset = NULL,
                         out = paste(file, ".subset", sep = "")) {
  
  sequences <- readDNAStringSet(file)
  if (all(as.character(subset) %in% names(sequences)) == FALSE)
    stop("There are names in 'subset' not present in the fasta file")
  pos <- match(as.character(subset), names(sequences))
  writeXStringSet(sequences[pos], filepath = out)
}


#' Parallelized NMDS Dimension Screening
#'
#' This function calculates and visualizes stress values for Non-Metric Multidimensional Scaling (NMDS)
#' across different numbers of dimensions (1 to 10). The stress calculation is parallelized to improve
#' performance, and the results are plotted to help choose the optimal number of dimensions.
#'
#' @param x A community data matrix (samples as rows, species as columns) or a distance matrix.
#' @param ncores An integer specifying the number of CPU cores to use for parallel processing.
#'   Defaults to `parallel::detectCores() - 1` (one less than the total number of available cores).
#'
#' @return A plot showing stress values for NMDS with dimensions ranging from 1 to 10. Each dimension
#'   is represented by 10 replicate stress values.
#'
#' @import vegan
#' @importFrom parallel mclapply detectCores
#' @importFrom graphics plot points
#' @export
#'
#' @examples
#' # Load example data
#' library(vegan)
#' data(varespec, package = "vegan")
#'
#' # Run NMDS dimension screening
#' nmds_screen_parallel(varespec, ncores = 2)

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

#' Standard Ordination Plot with ggplot2
#'
#' This function creates a standardized ordination plot using `ggplot2`.
#' It includes points, ellipses, and reference lines for visualizing NMDS or PCoA results. The function allows
#' customization of color, shape, and handling of missing values.
#'
#' @param .data A data frame containing ordination coordinates and additional metadata.
#' @param ordi A character string specifying the ordination method. Options are "NMDS" and "PcoA".
#' @param .color A column in `.data` to use for coloring points and ellipses. This should be a categorical variable.
#' @param .drop_na A column in `.data` used to filter out rows with missing values. Rows with `NA` in this column will be removed.
#'
#' @return A `ggplot` object representing the NMDS plot with points, ellipses, and reference lines.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import rlang
#' @export
#'
#' @examples
#' # Load example data
#' library(vegan)
#' data(dune, package = "vegan")
#' 
#' # Run NMDS
#' nmds_result <- metaMDS(dune, k = 2)
#' nmds_scores <- scores(nmds_result, display = "sites")
#' 
#' # Adding factor to data
#' dune$Management <- sample(c("A", "B"), nrow(dune), replace = TRUE)
#' nmds_data <- cbind(nmds_scores, dune)
#'
#' # Create NMDS plot
#' gg_ordi(nmds_data, .color = Management, ordi = "NMDS", .drop_na = Management)
# Standard Ordination plot
gg_ordi <- function(.data, .color, ordi, .drop_na = NULL) {
  # Input validation
  if (!inherits(.data, "data.frame")) {
    cli::cli_abort("{.arg .data} must be a data frame")
  }
  
  # Capture quosures safely
  .color <- rlang::enquo(.color)
  .drop_na <- rlang::enquo(.drop_na)
  
  ORDIS <- c("NMDS", "PCoA")
  
  ordi <- match.arg(ordi, ORDIS)
  
  # Base ggplot elements
  base_plot <- .data %>%
    {
      if (!rlang::quo_is_null(.drop_na)) {
        tidyr::drop_na(., !!.drop_na)
      } else {
        .
      }
    } %>%
    ggplot() +
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
    )
  
  # Add coordinates based on ordination type
  if (ordi == "NMDS") {
    base_plot <- base_plot +
      aes(x = NMDS1, y = NMDS2)
  } else if (ordi == "PCoA") {
    base_plot <- base_plot +
      aes(x = Dim1, y = Dim2)
  }
  
  # Add points and ellipses
  final_plot <- base_plot +
    geom_point(
      aes(color = !!.color),
      stroke = 1,
      alpha = 0.5,
      na.rm = TRUE
    ) +
    stat_ellipse(
      aes(color = !!.color),
      geom = "path",
      linewidth = 1.3,
      type = "t",
      level = 0.95,
      show.legend = TRUE
    )
  
  return(final_plot)
}
