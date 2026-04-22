#' Find Core Microbial Taxa
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
#' @importFrom dplyr sym left_join group_by summarise transmute arrange desc mutate n last select slice_head
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom rlang ensym as_name .data
#' @importFrom vegan decostand
#' @importFrom cli cli_text cli_warn cli_abort cli_alert_success cli_alert_info
#' @importFrom utils combn tail
#' @importFrom tidyr gather 
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
find_core <- function(ps_obj, 
                      priority_var = "sampling_date", 
                      nReads = 1000) {
    
    # 1. Data Extraction
    # Ensure we treat this as a matrix and keep row/col names
    otu <- as.matrix(as.data.frame(otu_table(ps_obj, taxa_are_rows = TRUE)))
    map <- as.data.frame(as.matrix(sample_data(ps_obj))) %>% 
        rownames_to_column("sequence_name")
    
    # 2. Occupancy and Ranking Logic
    otu_PA <- 1*((otu > 0) == 1)
    otu_occ <- rowSums(otu_PA) / ncol(otu_PA)
    otu_rel <- apply(vegan::decostand(otu, method = "total", MARGIN = 2), 1, mean)
    
    occ_abun <- data.frame(otu_occ = otu_occ, otu_rel = otu_rel) %>%
        rownames_to_column('otu')
    
    PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
        gather(sequence_name, abun, -otu) %>%
        left_join(map, by = 'sequence_name') %>%
        group_by(otu, !!sym(priority_var)) %>%
        summarise(time_freq = sum(abun > 0) / length(abun),
                  coreTime = ifelse(time_freq == 1, 1, 0), .groups = "drop") %>%
        group_by(otu) %>%
        summarise(sumF = sum(time_freq),
                  sumG = sum(coreTime),
                  nS = length(unique(!!sym(priority_var))),
                  Index = (sumF + sumG) / nS)
    
    # rank is the weighting index, represents how consistently an OTU appears 
    # across your priority_var
    otu_ranked <- occ_abun %>%
        left_join(PresenceSum, by = 'otu') %>%
        transmute(otu = otu, rank = Index) %>%
        arrange(desc(rank))
    
    # 3. Bray-Curtis Contribution (The Loop)
    BCaddition <- NULL
    
    # Crucial: Pre-calculate the combinations once
    # This uses the number of columns (samples) in your OTU table
    sample_combos <- utils::combn(ncol(otu), 2)
    x_names <- apply(sample_combos, 2, function(x) paste(colnames(otu)[x], collapse = ' - '))
    
    # Initialize with the 1st ranked OTU
    # drop = FALSE keeps it as a matrix so ncol() and combinations work
    start_matrix <- otu[otu_ranked$otu[1], , drop = FALSE]
    
    x <- apply(sample_combos, 2, function(idx) sum(abs(start_matrix[, idx[1]] - start_matrix[, idx[2]])) / (2 * nReads))
    BCaddition <- data.frame(x_names = x_names, "1" = x, check.names = FALSE)
    
    # Loop through remaining OTUs
    # Note: running this for ALL OTUs (nrow) can be very slow if you have thousands
    for(i in 2:nrow(otu_ranked)) {                   
        otu_add <- otu_ranked$otu[i]                        
        add_matrix <- otu[otu_add, , drop = FALSE]
        start_matrix <- rbind(start_matrix, add_matrix)
        
        # Calculate BC on the growing matrix
        x <- apply(sample_combos, 2, function(idx) sum(abs(start_matrix[, idx[1]] - start_matrix[, idx[2]])) / (2 * nReads))
        df_a <- data.frame(x_names = x_names, BC_val = x)
        names(df_a)[2] <- as.character(i)
        
        BCaddition <- left_join(BCaddition, df_a, by = 'x_names')
    }
    
    # 4. Process Results
    temp_BC_matrix <- BCaddition %>% 
        column_to_rownames("x_names") %>% 
        as.matrix()
    
    BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>% 
        gather(comparison, BC, -rank) %>%
        group_by(rank) %>%
        summarise(MeanBC = mean(BC)) %>%
        mutate(rank_num = as.numeric(as.character(rank))) %>%
        arrange(rank_num) %>%
        mutate(proportionBC = MeanBC / max(MeanBC))
    
    Increase <- BC_ranked$MeanBC[-1] / BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
    BC_ranked$IncreaseBC <- c(0, Increase)
    
    # 5. Elbow Method
    fo_difference <- function(pos, df) {
        if(pos == 1 || pos == nrow(df)) return(0)
        left <- (df$MeanBC[pos] - df$MeanBC[1]) / pos
        right <- (df$MeanBC[nrow(df)] - df$MeanBC[pos]) / (nrow(df) - pos)
        return(left - right)
    }
    BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference, df = BC_ranked)
    
    elbow_val <- which.max(BC_ranked$fo_diffs)
    lastCall_val <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])))
    
    # This is added by me to asses the core
    BC_ranked <- BC_ranked %>%
        mutate(
            rank_num = as.numeric(as.character(rank)),
            # Extract the 'otu' column from the row matching the rank_num
            otu_added = otu_ranked$otu[rank_num],
            # Calculate the % increase correctly
            delta_pct_max_BC = (IncreaseBC - 1) * 100 
        ) %>%
        # Relocate for better readability
        dplyr::select(rank, rank_num, otu_added, everything())
    
    # 6. Output as List
    return(list(
        otu_ranked = otu_ranked,
        BC_ranked_summary = BC_ranked,
        elbow_threshold = elbow_val,
        two_percent_threshold = lastCall_val,
        full_matrix = BCaddition
    ))
}
