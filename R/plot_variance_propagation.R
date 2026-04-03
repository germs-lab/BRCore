#' Variance propagation diagnostic for rarefaction
#'
#' This function evaluate the variance generated during multiple rarefaction iterations. 
#' It compares raw vs rarefied diversity metrics calculated at each iterations.
#'
#' @param physeq_obj Raw phyloseq object
#' @param rarefied Output from multi_rarefy(). Either a list of dataframes or and array. 
#' @param q Hill number order (q = 0 for richness, q = 1 for Shannon, q = 2 for Simpson)
#' @param distance "bray" or "jaccard" (beta only)
#' @param group_var A grouping variable to use gor grouping as in the sample_data()
#' @param group_color A color variable to use present in the sample_data()
#'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom phyloseq otu_table sample_data taxa_are_rows
#' @importFrom cli cli_h1 cli_h2 cli_alert_info cli_alert_warning cli_alert_success cli_alert_danger
#' @importFrom utils head
#' @importFrom vegan rrarefy
#' @importFrom dplyr left_join
#'  
#' @return ggplot object comparing raw vs rarefied diversity distributions across iterations.
#' 
#' @examples
#' \donttest{
#' library(phyloseq)
#' library(BRCore)
#' # Example using the bcse dataset, comparing hill q=1 between Poplar and Switchgrass plots
#' bcse_filt <- bcse %>% 
#' subset_samples(Crop %in% c("Poplar", "Switchgrass"))
#' bcse_rarefied_otutable_filt <-
#'  multi_rarefy(
#'        physeq = bcse_filt,
#'        depth_level = 1000,
#'        num_iter = 100,
#'        .as_array = FALSE,
#'        set_seed = 7643
#'    )
#'
#'plot_variance_propagation(
#'    physeq   = bcse_filt,
#'    rarefied = bcse_rarefied_otutable_filt,
#'    q        = 1,
#'    group_var = "Crop",
#'    group_color = "Plot"
#') + scale_color_viridis_d(option = "turbo")
#' }
#'
#' @export
plot_variance_propagation <- function(
        physeq_obj,
        rarefied,
        q = 0,
        group_var,
        group_color
) {
    cli::cli_h1("Rarefaction Variance Propagation Validation")

    # Extract raw OTU
    if (!inherits(physeq_obj, "phyloseq")) {
        cli::cli_abort(
            "{.arg physeq_obj} must be a 'phyloseq' object.\nYou've supplied a {class(physeq_obj)[1]} vector."
        )
    }
    
    cli::cli_alert_success("Input phyloseq object is valid!")
    
    otu_raw <- as(phyloseq::otu_table(physeq_obj), "matrix")
    if (phyloseq::taxa_are_rows(physeq_obj)) {
        otu_raw <- t(otu_raw)
    }
    
    meta <- data.frame(phyloseq::sample_data(physeq_obj))
    meta$sample_id <- rownames(meta)
    
    # Match samples to rarefied
    get_sample_ids <- function(x) {
        if (is.array(x)) rownames(x[, , 1])
        else if (is.list(x)) rownames(x[[1]])
        else rownames(x)
    }
    
    rare_samples <- get_sample_ids(rarefied)
    common_samples <- intersect(rownames(otu_raw), rare_samples)
    
    otu_raw <- otu_raw[common_samples, , drop = FALSE]
    meta    <- meta[common_samples, , drop = FALSE]
    
    cli::cli_alert_info("Hill number order selected, q= {.val {q}}")

    raw_vals <- apply(otu_raw, 1, .hill_number, q = q)
    
    raw_df <- data.frame(
        sample_id = names(raw_vals),
        value = raw_vals,
        method = "Raw"
    )
    
    # number of iterations in rarefied object
    get_n_iter <- function(x) {
        if (is.array(x)) return(dim(x)[3])
        if (is.list(x)) return(length(x))
        return(1)
    }
    
    n_iter <- get_n_iter(rarefied)
    
    cli::cli_alert_info("Number of rarefaction iterations, n_iter= {.val {n_iter}}")
    
    # replicate raw values
    raw_df <- raw_df[rep(seq_len(nrow(raw_df)), n_iter), ]
    
    # Rarefied alpha (ALL iterations)
    get_iter_list <- function(x) {
        if (is.array(x)) {
            lapply(seq_len(dim(x)[3]), function(i) x[, , i])
        } else if (is.list(x)) {
            x
        } else {
            list(x)
        }
    }
    
    iter_list <- get_iter_list(rarefied)
    
    rare_df <- lapply(seq_along(iter_list), function(i) {
        mat <- as.matrix(iter_list[[i]])
        
        vals <- apply(mat, 1, .hill_number, q = q)
        
        data.frame(
            sample_id = names(vals),
            value = vals,
            method = "Rarefied"
        )
    }) %>% bind_rows()
    
    # Combine
    all_df <- bind_rows(raw_df, rare_df) %>%
        left_join(meta, by = "sample_id")
    
    cli::cli_alert_info("Comparison plot generated!")
    
    p <- ggplot(all_df,
                aes(x = .data[[group_var]], y = value, color = .data[[group_color]])) +
        geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 0.8) +
        facet_wrap(~method) +
        theme_classic() +
        theme(
            strip.text = element_text(size = 10, face = "bold"),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 9),
            legend.background = element_blank(),
            legend.key.height = unit(0.4, "cm"),
            legend.key.width = unit(0.4, "cm")) +
        guides(color = guide_legend(
            override.aes = list(size = 5, shape = 15)  
        )) +
        labs(
            x = group_var,
            y = paste0("Hill number (q = ", q, ")"),
            title = "Raw vs Rarefied diversity"
        )
    
    return(p)
}

#' Helper function to calculate hill numbers 
#' @param x Sample seqeunce counts
#' @param q Hill number order (q = 0 for richness, q = 1 for Shannon, q = 2 for Simpson)
.hill_number <- function(x, q) {
    p <- x / sum(x)
    p <- p[p > 0]
    
    if (q == 0) return(length(p))
    if (q == 1) return(exp(-sum(p * log(p))))
    if (q == 2) return(1 / sum(p^2))
    
    (sum(p^q))^(1 / (1 - q))
}
