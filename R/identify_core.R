#' Identify Core Microbial Taxa
#'
#' @description
#' This function identifies core microbial taxa based on abundance-occupancy
#' distributions and their contributions to Bray-Curtis similarity between
#' biological samples. Core taxa are selected using either a "last % increase"
#' or "elbow" method implementing the method developed by Shade
#' and Stopnisek (2019) Curr Opin Microbiol, see below for details.
#'
#' @param physeq_obj A `phyloseq` object with at least `otu_table` and
#' `sample_data`.
#' @param priority_var The column name in the `sample_data` (e.g.
#' sampling_date", "site") that is used for prioritizing the core microbiome.
#' @param increase_value Increase value (numeric, scalar) used in the
#' calculation (default 0.02) for "increase". The "elbow" is always calculated and returned as \code{elbow_core} (see below for details).
#' @param abundance_weight Numeric in `[0,1]`; how much to weight mean relative
#' abundance in the ranking score. `0` (default) uses occupancy/composite only.
#' `1` ranks purely by abundance. Values in between blend the two (e.g.,
#' abundance_weight = 0.3 gives 70% occupancy/composite + 30% abundance).
#' @param max_otus Optional integer to limit analysis to the top N ranked OTUs.
#' If NULL (default), all OTUs are analyzed. Useful for large datasets
#' (>5000 OTUs)
#' @param depth_level Integer. The sequencing depth used for normalization in
#' Bray-Curtis calculations. If data is rarefied, this is automatically set
#' to the rarefaction depth. For unrarefied data, samples with depth below
#' this threshold are excluded from pairwise comparisons.
#' @param num_iter Integer. Number of subsampling iterations used when
#' calculating average dissimilarity for unrarefied data (default 100).
#' Ignored if data is already rarefied.
#' @param seed Optional integer to set the RNG seed for reproducibility.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{bray_curtis_ranked} tibble with `rank`, mean Bray-Curtis
#'  similarity across sample pairs (`MeanBC`) at each cumulative `rank`,
#'  normalized proportion (`proportionBC`), the multiplicative `IncreaseBC`,
#' and the elbow metric (`elbow_slope_diffs`).
#' (`proportionBC`), the multiplicative `IncreaseBC`, and the elbow metric
#' (`elbow_slope_diffs`).
#'   \item \code{otu_ranked} tibble with ranked OTU/ASVs .
#'   \item \code{abundance_occupancy} tibble with OTU/ASVs names, occupancy
#'      (`otu_occ`), and mean relative abundance (`otu_rel`).
#'   \item \code{priority_var} character, the variable used for prioritizing
#' the core.
#'   \item \code{elbow} core set identified by elbow method (integer).
#'   \item \code{bc_increase} core set identified by last % BC-increase
#' (integer).
#'   \item \code{increase_value} increase value (numeric, scalar) used in the
#' calculation (e.g. 0.02).
#'   \item \code{elbow_core} core OTU/ASVs using elbow method (character
#' vector).
#'   \item \code{increase_core} core OTU/ASVs using last % BC-increase method
#'  (character vector).
#'   \item \code{otu_table} otu_table counts (otu x samples) used (data.frame).
#'   \item \code{sample_metadata} samples metadata (data.frame).
#'   \item \code{taxonomy_table} taxonomy if present (data.frame); otherwise
#'  NULL.
#' }
#'
#' @details
#' The core set is defined using two separate methods:
#'
#' The function rank OTU/ASVs by occupancy (optionally with abundance
#' weighting: `rank_score = (1 - weight) * Index + weight * scaled_abundance`,
#' where scaled_abundance is mean relative abundance rescaled to `[0,1]`). For
#' each `k = 1..K`, recompute `S_k` as the mean Bray-Curtis similarity across
#' all sample pairs using only the first `k` ranked OTUs; when `k = K`, this
#' yields `S_K`, the value computed with all OTUs. Normalizing by `S_K`
#' gives `C_k = S_k / S_K`.
#'
#' The **elbow** is the point of diminishing returns: for each `k`, compare the
#' average *left* slope `(S_k - S_1) / (k - 1)` to the average *right* slope
#' `(S_K - S_k) / (K - k)`, and choose the `k` that maximizes `(left - right)`.
#'
#' The **last percent Bray-Curtis increase** method uses the same accumulation
#' curve, examine the multiplicative step when adding the `k-th` OTU:
#' `Increase_k = S_k / S_{k-1}` (equivalently, `Increase_k = C_k / C_{k-1}`).
#' Choose the largest `k` such that `Increase_k >= 1 + p`, where `p` is your
#' chosen percent threshold (increase_value; recommended `p >= 0.02` or `2%`).
#' This selects the final rank for which adding one more taxon still increases
#' the explained Bray-Curtis similarity by at least `p`.
#'
#' @references Shade A, Stopnisek N (2019) Abundance-occupancy
#' distributions to prioritize plant core microbiome membership. Current
#' Opinion in Microbiology, 49:50-58
#' doi:https://doi.org/10.1016/j.mib.2019.09.008
#'
#' @section Dependencies:
#' Requires \pkg{phyloseq}, \pkg{dplyr}, \pkg{tidyr}, \pkg{tibble}, \pkg{rlang},
#' and \pkg{vegan}.
#'
#' @importFrom phyloseq sample_sums taxa_are_rows otu_table sample_data tax_table
#' @importFrom dplyr left_join group_by summarise transmute arrange desc mutate n last select slice_head
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom rlang ensym as_name .data
#' @importFrom vegan decostand
#' @importFrom cli cli_text cli_warn cli_abort cli_alert_success cli_alert_info
#' @importFrom utils combn tail
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' library(BRCore)
#' # Example using your switchgrass phyloseq object and grouping variable
#' # 'sampling_date'
#' data("switchgrass", package = "BRCore")
#'
#' res <- identify_core(
#'   physeq_obj     = switchgrass,
#'   priority_var   = "sampling_date",
#'   increase_value = 0.02,
#'   seed           = 091825
#' )
#'
#' # Inspect results
#' str(res)
#' }
#'
#' @export
identify_core <- function(
  physeq_obj,
  priority_var,
  increase_value = 0.02,
  abundance_weight = 0,
  max_otus = NULL,
  depth_level = 1000,
  num_iter = 100,
  seed = NULL
) {
  # input checks ----
  cli::cli_text("\nSeed used: {seed}\n")

  if (is.null(seed)) {
    cli::cli_warn("No seed was set. Results may not be reproducible.")
  } else {
    set.seed(seed)
  }

  .phyloseq_class_check(physeq_obj)

  # define arguments ----

  # Check if samples are rarefied (all have same depth, accounting for floating-point precision)
  min_sum <- min(sample_sums(physeq_obj))
  max_sum <- max(sample_sums(physeq_obj))
  is_rarefied <- abs(min_sum - max_sum) < 1e-6

  if (is_rarefied) {
    depth_level <- round(min_sum) # Use rounded value for display
    cli::cli_alert_info(
      "otu_table() is rarefied at a depth of: {.val {depth_level}}"
    )
  } else {
    cli::cli_alert_warning(
      "The otu_table() is not rarefied! \n Using depth_level={.val {depth_level}} for rarefaction and normalization in Bray-Curtis calculations. \n Adjust depth_level according to your objectives."
    )
  }

  otu <- .extract_otu_matrix(physeq_obj, samples_as_rows = FALSE)
  map <- sample_data(physeq_obj) |> as("data.frame")
  map$sample_id <- rownames(map)

  ### check for a tax_table if present, if not just print a warning and continue.
  taxa <- NULL
  tx <- tryCatch(tax_table(physeq_obj), error = function(e) NULL)

  if (!is.null(tx)) {
    txm <- tryCatch(as(tx, "matrix"), error = function(e) NULL)
    if (!is.null(txm) && is.matrix(txm) && ncol(txm) > 0L) {
      taxa <- as.data.frame(txm, stringsAsFactors = FALSE)
    }
  }

  if (is.null(taxa)) {
    cli::cli_alert_info(
      "No taxonomy found (or empty). Continuing without taxonomy."
    )
  }

  # core prioritizing variable ----
  data_var <- priority_var
  cli::cli_alert_success("Core prioritizing variable: {data_var}")

  # validate abundance_weight ----
  if (
    !is.numeric(abundance_weight) ||
      length(abundance_weight) != 1L ||
      is.na(abundance_weight)
  ) {
    stop("`abundance_weight` must be a single numeric in [0,1].")
  }
  if (abundance_weight < 0 || abundance_weight > 1) {
    cli::cli_warn(
      "`abundance_weight`={abundance_weight} is outside [0,1]; clamping."
    )
    abundance_weight <- max(0, min(1, abundance_weight))
  }

  # abundance occupancy ----
  # occupancy and mean rel. abundance
  otu_PA <- 1 * (otu > 0)
  otu_occ <- rowSums(otu_PA) / ncol(otu_PA)
  otu_rel <- apply(
    decostand(otu, method = "total", MARGIN = 2),
    1,
    mean
  )

  occ_abun <- data.frame(otu_occ = otu_occ, otu_rel = otu_rel) |>
    tibble::rownames_to_column("otu")

  # PresenceSum ----

  PresenceSum <- data.frame(
    otu = as.factor(rownames(otu)),
    otu,
    check.names = FALSE
  ) |>
    pivot_longer(-otu, names_to = "sample_id", values_to = "abun") |>
    left_join(map, by = "sample_id") |>
    group_by(otu, .data[[data_var]]) |>
    summarise(
      time_freq = sum(abun > 0) / length(abun), # frequency within date
      coreTime = ifelse(time_freq == 1, 1, 0), # 1 if occupancy=1 at that date
      .groups = "drop"
    ) |>
    group_by(otu) |>
    summarise(
      sumF = sum(time_freq),
      sumG = sum(coreTime),
      nS = 2L * length(.data[[data_var]]), # K per original
      Index = (sumF + sumG) / nS, # can be up to 2 (original behavior)
      .groups = "drop"
    )

  # Rank and weight in abundance ----
  otu_ranked <- occ_abun |>
    left_join(PresenceSum, by = "otu")

  rank_method <- if (abundance_weight > 0) {
      occupancy_weight <- 1 - abundance_weight
      paste0("Occupancy: ", occupancy_weight, " and Abundance: ", abundance_weight) 
  } else {
      "Index only"
  }

  if (abundance_weight == 0) {
    # No abundance weighting - use Index directly
    otu_ranked <- otu_ranked |>
      mutate(rank = Index) |>
      arrange(desc(rank))

    cli::cli_alert_info(
      "Ranked by {rank_method}"
    )
  } else {
    otu_ranked <- otu_ranked |>
      mutate(
        occ_norm = otu_occ / max(otu_occ, na.rm = TRUE),
        abun_norm = otu_rel / max(otu_rel, na.rm = TRUE),
        spatial_weight = (abundance_weight * abun_norm) +
          ((1 - abundance_weight) * occ_norm),
        rank = spatial_weight * Index
      ) |>
      arrange(desc(rank), desc(Index), otu) |>
      select(
        otu,
        rank,
        Index,
        spatial_weight,
        occ_norm,
        abun_norm,
        otu_occ,
        otu_rel
      )
    cli::cli_alert_info(
      "Ranked by {rank_method}"
    )
  }

  # Optional: Limit analysis to top N OTUs
  if (!is.null(max_otus)) {
    if (!is.numeric(max_otus) || max_otus < 1) {
      cli::cli_abort("`max_otus` must be a positive integer.")
    }

    max_otus <- as.integer(max_otus)

    if (max_otus < nrow(otu_ranked)) {
      total_otus <- nrow(otu_ranked)
      otu_ranked <- otu_ranked |>
        slice_head(n = max_otus)

      cli::cli_alert_info(
        "Limiting analysis to top {max_otus} of {total_otus} OTUs (ranked by {rank_method})."
      )
    } else {
      cli::cli_alert_info(
        "max_otus ({max_otus}) >= total OTUs ({nrow(otu_ranked)}). Using all OTUs (ranked by {rank_method})."
      )
    }
  }

  otu_ranked_ordered <- otu_ranked$otu

  # BC accumulation ----
  # cumulative BC across samples while adding taxa in rank order
  # pairwise BC on *current* subset of taxa, normalized by depth_level, matching your formula

  # start with the first ranked OTU
  cli::cli_alert_info(
    "Ranking OTUs based on BC dissimilarity, starting at {Sys.time()}"
  )

  start_idx <- match(otu_ranked$otu[1], rownames(otu))
  start_matrix <- matrix(
    otu[start_idx, ],
    nrow = 1,
    dimnames = list(otu_ranked$otu[1], colnames(otu))
  )

  bc_vec <- .calculate_bc(
    start_matrix,
    depth_level,
    num_iterations = num_iter,
    is_rarefied = is_rarefied
  )

  BCaddition <- data.frame(
    x_names = bc_vec$names,
    `1` = bc_vec$values,
    check.names = FALSE
  )

  if (nrow(otu_ranked) > 1) {
    progressbar_calc_bc <- cli::cli_progress_bar(
      name = "Calculating BC rankings",
      total = nrow(otu_ranked) - 1,
      format = "{cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
      .auto_close = TRUE
    )

    for (i in 2:nrow(otu_ranked)) {
      add_idx <- match(otu_ranked$otu[i], rownames(otu))
      add_matrix <- matrix(
        otu[add_idx, ],
        nrow = 1,
        dimnames = list(otu_ranked$otu[i], colnames(otu))
      )
      start_matrix <- rbind(start_matrix, add_matrix)

      bc_vec <- .calculate_bc(
        start_matrix,
        depth_level,
        num_iterations = num_iter,
        is_rarefied = is_rarefied
      )

      BCaddition <- left_join(
        BCaddition,
        data.frame(
          x_names = bc_vec$names,
          value = bc_vec$values,
          check.names = FALSE
        ),
        by = "x_names"
      )
      names(BCaddition)[ncol(BCaddition)] <- as.character(i)

      cli::cli_progress_update()
    }
  }

  temp_BC_matrix <- BCaddition |>
    tibble::column_to_rownames("x_names") |>
    as.matrix()

  BC_ranked <- data.frame(
    rank = as.factor(rownames(t(temp_BC_matrix))),
    t(temp_BC_matrix),
    check.names = FALSE
  ) |>
    pivot_longer(
      -rank,
      names_to = "comparison",
      values_to = "BC"
    ) |>
    group_by(.data$rank) |>
    summarise(MeanBC = mean(.data$BC), .groups = "drop") |>
    arrange(MeanBC) |>
    mutate(proportionBC = .data$MeanBC / max(.data$MeanBC)) |>
    mutate(
        pct_change_proportionBC = 100 * (proportionBC / lag(proportionBC) - 1)
        )

  # increase method: multiplicative increase between successive ranks ----
  if (nrow(BC_ranked) >= 2) {
    Increase <- BC_ranked$MeanBC[-1] / BC_ranked$MeanBC[-nrow(BC_ranked)]
    increaseDF <- data.frame(
      IncreaseBC = c(0, Increase),
      rank = factor(seq_len(length(Increase) + 1))
    )
  } else {
    increaseDF <- data.frame(IncreaseBC = 0, rank = factor(1))
  }
  BC_ranked <- left_join(BC_ranked, increaseDF, by = "rank")

  # elbow method: by forward-backward slope difference ----
  elbow_slope_differences <- function(pos) {
    left <- (BC_ranked$MeanBC[pos] - BC_ranked$MeanBC[1]) / pos
    right <- (BC_ranked$MeanBC[nrow(BC_ranked)] - BC_ranked$MeanBC[pos]) /
      max(1, (nrow(BC_ranked) - pos))
    left - right
  }
  BC_ranked$elbow_slope_diffs <- vapply(
    seq_len(nrow(BC_ranked)),
    elbow_slope_differences,
    numeric(1)
  )

  elbow <- which.max(BC_ranked$elbow_slope_diffs)

  # FALLBACK: if elbow is empty, NA, or 0, default to 1
  if (length(elbow) == 0 || is.na(elbow) || elbow == 0) {
    cli::cli_warn(
      "Elbow calculation failed (possibly due to low variation). Defaulting to 1."
    )
    elbow <- 1
  }

  # threshold cut based on multiplicative increase
  thr <- 1 + increase_value
  valid_increases <- which(BC_ranked$IncreaseBC >= thr)

  lastCall <- if (length(valid_increases)) {
    utils::tail(valid_increases, 1)
  } else {
    1
  }

  # identified core sets ----

  elbow_core <- otu_ranked_ordered[seq_len(elbow)]
  cli::cli_alert_success(
    "Elbow method identified {.val {length(elbow_core)}} core OTUs"
  )
  increase_core <- otu_ranked_ordered[seq_len(lastCall)]
  cli::cli_alert_success(
    "% increase method identified {.val {length(increase_core)}} core OTUs"
  )

  # return ----
  out <- list(
    bray_curtis_ranked = BC_ranked,
    otu_ranked = otu_ranked,
    abundance_occupancy = occ_abun,
    priority_var = data_var,
    abundance_weight = abundance_weight,
    elbow = as.integer(elbow),
    bc_increase = as.integer(lastCall),
    increase_value = increase_value,
    elbow_core = elbow_core,
    increase_core = increase_core,
    otu_table = otu,
    metadata = map,
    taxonomy = taxa
  )

  class(out) <- c("identify_core_result", class(out))

  cli::cli_alert_success("Analysis complete!")

  out
}


#' Calculate Bray-Curtis Dissimilarity Between Sample Pairs
#'
#' This function calculates the Bray-Curtis dissimilarity between all pairs of samples in a matrix.
#' The dissimilarity is normalized by the total number of reads (`depth_level`) to account for differences
#' in sequencing depth.
#'
#' @param matrix A numeric matrix or data frame where rows represent taxa (e.g., ASVs) and columns
#'   represent samples. The column names of the matrix are used to generate pairwise sample names.
#' @param depth_level A numeric value representing the total number of reads used for normalization.
#' @param num_iterations
#'   Typically, this is the minimum or average sequencing depth across samples.
#' @return A list containing:
#'   \itemize{
#'     \item \code{values}: A numeric vector of Bray-Curtis dissimilarity values for each sample pair.
#'     \item \code{names}: A character vector of pairwise sample names (e.g., "Sample1-Sample2").
#'   }
#' @export
#'
#' @noRd
#' @keywords internal

.calculate_bc <- function(
  matrix,
  depth_level,
  num_iterations,
  is_rarefied = TRUE
) {
  if (nrow(matrix) == 0) {
    cli::cli_alert_warning("{.arg matrix} is empty. Enter a non-empty matrix.")
    return(list(values = numeric(0), names = character(0)))
  }

  if (ncol(matrix) < 2) {
    cli::cli_alert_warning(
      "{.arg matrix} has fewer than 2 columns. Need at least 2 columns to calculate pairwise distances."
    )
    return(list(values = numeric(0), names = character(0)))
  }

  # Conditional handling for rarefied vs unrarefied data
  if (is_rarefied) {
    # Data already rarefied - use vegdist directly without subsampling
    sample_pairs <- utils::combn(ncol(matrix), 2)
    pair_labels <- apply(sample_pairs, 2, function(x) {
      paste(colnames(matrix)[x], collapse = " - ")
    })

    bc_vegan <- as.vector(vegan::vegdist(t(matrix), method = "bray"))

    # Per-pair normalization vegdist used
    pair_sums <- apply(sample_pairs, 2, function(x) {
      sum(matrix[, x[1]]) + sum(matrix[, x[2]])
    })

    # Reverse vegdist normalization and apply our own
    bc_values <- bc_vegan * pair_sums / (2 * depth_level)
  } else {
    # Unrarefied data - filter samples by depth FIRST
    floating_depth <- depth_level * 0.9999999 # account for floating-point precision

    keep_samples <- colSums(matrix) >= floating_depth #ASV/OTUs should be columns

    filtered_matrix <- matrix[, keep_samples, drop = FALSE]

    if (ncol(filtered_matrix) < 2) {
      return(list(values = numeric(0), names = character(0)))
    }

    # NOW compute sample pairs from the filtered matrix
    sample_pairs <- utils::combn(ncol(filtered_matrix), 2)
    pair_labels <- apply(sample_pairs, 2, function(x) {
      paste(colnames(filtered_matrix)[x], collapse = " - ")
    })

    # Use avgdist with rarefaction
    bc_values <- as.vector(
      suppressWarnings(
        vegan::avgdist(
          t(filtered_matrix),
          sample = depth_level,
          iterations = num_iterations,
          dmethod = "bray"
        )
      )
    )
  }

  list(
    values = bc_values,
    names = pair_labels
  )
}
