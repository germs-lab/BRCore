#' Calculate and append pre-rarefaction statistics to microbiome data
#'
#' This function adds read count, singleton count, Good's coverage,
#' and marks outlier samples to a `phyloseq` object or `data.frame`
#' based on the OTU/ASV abundance table.
#'
#' @param data A `phyloseq` object or a `data.frame` with samples
#'   as rows and taxa as columns.
#'
#' @return The same object (`phyloseq` or `data.frame`) with new
#'   columns:
#' \itemize{
#'   \item `read_num`
#'   \item `singleton_num`
#'   \item `goods_cov`
#'   \item `outlier`
#' }
#'
#' @details About Good's coverage. Initially developed by Alan
#'   Turing and I.J. Good during their cryptographic analyses in
#'   World War II, it was later adopted by ecologists, particularly
#'   in microbial diversity studies, to assess the completeness of
#'   a sample's representation of the overall community.
#' It's calculated as `1 - (F1/N)`, where `F1` is the number of
#'   OTUs (Operational Taxonomic Units) represented by only one
#'   individual (singletons) and `N` is the total number of
#'   individuals in the sample. For example, a Good's coverage of
#'   0.95, means that 5% of the reads in that sample are from OTUs
#'   that appear only once.
#'
#' @seealso [plot_rarefaction_metrics()] for visualizing these metrics, and
#' [multi_rarefy()] for performing rarefaction on a `phyloseq` object.
#'
#' @examples
#' library(phyloseq)
#' library(BRCore)
#'
#' data("bcse", package = "BRCore")
#'
#' # Adding metrics to a "phyloseq" object
#' bcse_metrics <- add_rarefaction_metrics(data = bcse)
#' sample_data(bcse_metrics)|>
#' head(10)
#'
#'
#' # Adding metrics to a "data.frame" count table object
#' bcse_otutable <- as.data.frame(
#'   as.matrix(otu_table(bcse))
#' )
#'
#' bcse_otutable_metrics <- add_rarefaction_metrics(
#'   data = bcse_otutable
#' )
#' bcse_otutable_metrics[
#'   head(seq_len(nrow(bcse_otutable_metrics)), 10),
#'   tail(seq_len(ncol(bcse_otutable_metrics)), 20)
#' ]
#'
#'
#' @importFrom phyloseq otu_table sample_data sample_data<- taxa_are_rows
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr group_by summarize mutate ungroup left_join
#' @importFrom tidyr pivot_longer
#' @importFrom stats IQR quantile
#' @importFrom cli cli_abort
#'
#' @export
add_rarefaction_metrics <- function(data) {
  # Internal function to identify outliers
  find_outlier <- function(x) {
    x < quantile(x, 0.25) - 1.5 * IQR(x) |
      x > quantile(x, 0.75) + 1.5 * IQR(x)
  }

  # Extract OTU table and transpose if needed
  if (inherits(data, "phyloseq")) {
    otu_mat <- .extract_otu_matrix(data, samples_as_rows = FALSE)
    otu_df <- as.data.frame(otu_mat)
  } else if (inherits(data, "data.frame")) {
    otu_df <- as.data.frame(t(data))
  } else {
    cli::cli_abort("Input must be a phyloseq object or a data.frame")
  }

  # Compute metrics
  df_stats <- otu_df |>
    rownames_to_column("otu_id") |>
    pivot_longer(
      -otu_id,
      names_to = "sample_id",
      values_to = "seq_num"
    ) |>
    group_by(.data$sample_id) |>
    summarize(
      read_num = sum(.data$seq_num),
      singleton_num = sum(.data$seq_num == 1),
      goods_cov = 100 * (1 - .data$singleton_num / .data$read_num)
    ) |>
    mutate(
      outlier = ifelse(
        find_outlier(log10(.data$read_num)),
        .data$read_num,
        NA
      )
    ) |>
    ungroup()

  # Append to sample metadata or return updated data.frame
  if (inherits(data, "phyloseq")) {
    sample_df <- sample_data(data) |>
      as.matrix() |>
      as.data.frame() |>
      rownames_to_column("sample_id")

    sample_df_updated <- left_join(
      sample_df,
      df_stats,
      by = "sample_id"
    ) |>
      as.data.frame() |>
      column_to_rownames("sample_id")

    sample_data(data) <- sample_df_updated
    return(data)
  } else {
    data_out <- data |>
      rownames_to_column("sample_id") |>
      left_join(df_stats, by = "sample_id") |>
      column_to_rownames("sample_id")

    return(data_out)
  }
}
