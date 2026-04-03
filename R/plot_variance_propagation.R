#' Variance propagation diagnostic for rarefaction
#'
#' Compare raw vs rarefied diversity distributions across iterations.
#'
#' @param physeq Raw phyloseq object
#' @param rarefied Output from iterated_rarefy_phyloseq()
#' @param measure "alpha" or "beta"
#' @param q Hill number order (alpha only)
#' @param distance "bray" or "jaccard" (beta only)
#' @param group_var grouping variable in sample_data
#'
#' @return ggplot
#' @export
plot_variance_propagation <- function(
        physeq,
        rarefied,
        q = 0,
        group_var,
        group_color
) {
    
    require(phyloseq)
    require(dplyr)
    require(ggplot2)
    
    # -----------------------------
    # Hill function
    # -----------------------------
    hill_number <- function(x, q) {
        p <- x / sum(x)
        p <- p[p > 0]
        
        if (q == 0) return(length(p))
        if (q == 1) return(exp(-sum(p * log(p))))
        if (q == 2) return(1 / sum(p^2))
        
        (sum(p^q))^(1 / (1 - q))
    }
    
    # -----------------------------
    # Extract raw OTU
    # -----------------------------
    otu_raw <- as(phyloseq::otu_table(physeq), "matrix")
    if (phyloseq::taxa_are_rows(physeq)) {
        otu_raw <- t(otu_raw)
    }
    
    meta <- data.frame(phyloseq::sample_data(physeq))
    meta$sample_id <- rownames(meta)
    
    # -----------------------------
    # Match samples to rarefied
    # -----------------------------
    get_sample_ids <- function(x) {
        if (is.array(x)) rownames(x[, , 1])
        else if (is.list(x)) rownames(x[[1]])
        else rownames(x)
    }
    
    rare_samples <- get_sample_ids(rarefied)
    common_samples <- intersect(rownames(otu_raw), rare_samples)
    
    otu_raw <- otu_raw[common_samples, , drop = FALSE]
    meta    <- meta[common_samples, , drop = FALSE]
    
    # -----------------------------
    # RAW alpha
    # -----------------------------
    raw_vals <- apply(otu_raw, 1, hill_number, q = q)
    
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
    
    # replicate raw values
    raw_df <- raw_df[rep(seq_len(nrow(raw_df)), n_iter), ]
    
    # -----------------------------
    # Rarefied alpha (ALL iterations)
    # -----------------------------
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
        
        vals <- apply(mat, 1, hill_number, q = q)
        
        data.frame(
            sample_id = names(vals),
            value = vals,
            method = "Rarefied"
        )
    }) %>% bind_rows()
    
    # -----------------------------
    # Combine
    # -----------------------------
    all_df <- bind_rows(raw_df, rare_df) %>%
        left_join(meta, by = "sample_id")
    
    # -----------------------------
    # Plot
    # -----------------------------
    p <- ggplot(all_df,
                aes(x = .data[[group_var]], y = value, color = .data[[group_color]])) +
        geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 0.8) +
        facet_wrap(~method) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
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

