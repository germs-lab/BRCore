#' Run multiple rarefaction for microbiome count tables
#'
#' This function performs multiple rarefaction on a `phyloseq` object by randomly
#' sub-sampling OTUs/ASVs within samples without replacement. The process is
#' repeated for a specified number of iterations, and the results are averaged.
#' Samples with fewer OTUs/ASVs than the specified `depth_level` are discarded.
#'
#' @param physeq A `phyloseq` object containing an OTU/ASV table.
#' @param depth_level An integer specifying the sequencing depth (number of
#'   OTUs/ASVs) to which samples should be rarefied.
#' @param num_iter An integer specifying the number of iterations to perform
#'   for rarefaction.
#' @param threads Number of threads (default = 4).
#' @param set_seed An optional integer to set the random seed for reproducibility (default = NULL).
#'
#' @return A data frame with taxa as rows and samples as columns. The values
#'   represent the average sequence counts calculated across all iterations.
#'   Samples with less than `depth_level` sequences are discarded.
#'
#'
#' @importFrom parallelly availableCores
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @importFrom dplyr  bind_rows group_by summarise across everything filter
#' @importFrom dplyr where
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom phyloseq otu_table
#' @importFrom vegan rrarefy
#' @importFrom cli cli_text cli_warn
#' @importFrom utils head
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' library(BRCore)
#' data("bcse", package = "BRCore")
#'
#' # Example rarefaction (single iteration, single core to keep examples fast)
#' otu_table_rare <- multi_rarefy(
#'    physeq = bcse,
#'    depth_level = 1000,
#'    num_iter = 100,
#'    threads = 2,
#'    set_seed = 7642
#')
#'
#' rowSums(otu_table_rare)
#' }
#'
#' @export
multi_rarefy <- function(
    physeq,
    depth_level,
    num_iter = 100,
    threads = get_available_cores(),
    set_seed = NULL
) {
    # Input validation ----
    cli::cli_text("\nSeed used: {set_seed}\n")

    if (is.null(set_seed)) {
        cli::cli_warn("No seed was set. Results may not be reproducible.")
    } else {
        set.seed(set_seed)
    }

    # Check object class
    if (!inherits(physeq, "phyloseq")) {
        stop("Input must be a phyloseq object, not a data.frame")
    }

    # Prepare data ----
    dataframe <- as.data.frame(
        as.matrix(t(otu_table(physeq, taxa_are_rows = TRUE)))
    )

    ### debug ###########################################

    printAndReturn <- function(x) {
        print("\n--- rowSums dplyr::across ---\n")
        print(x)
        print("\n--- before rounding ---\n")
        print(depth_level)
        print(class(x))
        print("\n--- x >= depth_level? ---\n")
        print(x >= depth_level)
        print("\n--- x formatted ---\n")
        print(format(x, nsmall = 20))
        print("\n--- x as integer?---\n")
        print(as.integer(x))
        print(as.integer(depth_level))
        print(as.integer(x) >= as.integer(depth_level))

        print("\n--- After rounding ---\n")
        y <- round(x)
        print(depth_level)
        print(class(y))
        print("\n--- y >= depth_level? ---\n")
        print(y >= depth_level)
        print("\n--- y formatted ---\n")
        print(format(y, nsmall = 20))
        print("\n--- y as integer?---\n")
        print(as.integer(y))
        print(as.integer(depth_level))
        print(as.integer(y) >= as.integer(depth_level))

        x
    }

    cat("taxa_are_rows(physeq):", phyloseq::taxa_are_rows(physeq), "\n")
    cat("otu_mat dim:", paste(dim(dataframe), collapse = " x "), "\n")
    cat(
        "nsamples:",
        phyloseq::nsamples(physeq),
        " ntaxa:",
        phyloseq::ntaxa(physeq),
        "\n"
    )

    # what you currently build:
    #dataframe <- as.data.frame(as.matrix(t(phyloseq::otu_table(physeq, taxa_are_rows = TRUE))))
    cat("dataframe dim:", paste(dim(dataframe), collapse = " x "), "\n")
    cat(
        "dataframe rownames head:",
        paste(head(rownames(dataframe)), collapse = ", "),
        "\n"
    )
    cat(
        "dataframe colnames head:",
        paste(head(colnames(dataframe)), collapse = ", "),
        "\n"
    )

    # sanity: do rows look like samples?
    cat(
        "Rows match sample_names?:",
        all(rownames(dataframe) %in% phyloseq::sample_names(physeq)),
        "\n"
    )
    cat(
        "Cols match taxa_names?:",
        all(colnames(dataframe) %in% phyloseq::taxa_names(physeq)),
        "\n"
    )

    # critical: row totals vs sample sums
    print(summary(rowSums(dataframe)))

    ### end debug #######################################

    # Parallel setup ----
    threads <- min(threads, availableCores())
    cl <- makeCluster(threads)
    on.exit(stopCluster(cl), add = TRUE)

    # Export needed objects/packages to workers
    clusterExport(
        cl,
        varlist = c("dataframe", "depth_level", ".single_rarefy", "set_seed"),
        envir = environment()
    )

    # Run rarefactions in parallel ----
    com_iter <- parLapply(cl, 1:num_iter, function(i) {
        set.seed(set_seed + i) # each worker gets a different but reproducible seed
        df <- .single_rarefy(dataframe, sample = depth_level)
        df <- as.data.frame(df)
        rownames_to_column(df, "sample_id")
    })

    # Aggregate results ----
    mean_data <- bind_rows(com_iter) |>
        group_by(sample_id) |>
        summarise(across(everything(), mean), .groups = "drop") |>
        filter(
            printAndReturn(round(rowSums(across(where(is.numeric))))) >=
                depth_level
        ) |>
        column_to_rownames("sample_id")

    # Remove ASVs/OTUs with zero total abundance
    mean_data <- mean_data[, colSums(mean_data) > 0]

    return(mean_data)
}


#' Single rarefaction iteration
#'
#' @param x A matrix with samples as rows and taxa as columns.
#' @param sample_size The number of sequences to sample.
#'
#' @return A matrix with rarefied counts.
#'
#' @noRd
#' @keywords internal
.single_rarefy <- function(x, sample_size) {
    x <- as.matrix(x)

    out <- t(apply(x, 1, function(row) {
        total <- sum(row)
        if (total < sample_size) {
            return(rep(NA_integer_, length(row)))
        }

        sampled <- sample(
            rep(seq_along(row), times = row),
            size = sample_size,
            replace = FALSE
        )
        tabulate(sampled, nbins = length(row))
    }))

    colnames(out) <- colnames(x)
    rownames(out) <- rownames(x)
    out
}
