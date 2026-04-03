#' Iterated rarefaction for phyloseq objects
#'
#' Performs one or multiple rarefactions using vegan::rrarefy().
#'
#' @param physeq A phyloseq object containing an otu_table.
#' @param depth Rarefaction depth.
#' @param n_iter Number of rarefaction iterations.
#' @param .array Logical; if TRUE return a 3D array, otherwise a list of data.frames.
#'
#' @return
#' - If n_iter = 1: a data.frame (samples x taxa)
#' - If n_iter > 1 and .array = TRUE: a 3D array (samples x taxa x iterations)
#' - If n_iter > 1 and .array = FALSE: a list of data.frames
#'
#' @importFrom vegan rrarefy
#' @export
multi_rarefy_new <- function(physeq, 
                             depth, 
                             n_iter = 1, 
                             .array = TRUE) {
    
    # --- checks ---
    if (!inherits(physeq, "phyloseq")) {
        stop("physeq must be a phyloseq object")
    }
    
    if (n_iter < 1) {
        stop("n_iter must be >= 1")
    }
    
    # --- extract OTU table safely ---
    otu <- phyloseq::otu_table(physeq)
    
    # convert to matrix FIRST
    otu_mat <- as(otu, "matrix")
    
    # enforce samples as rows
    if (phyloseq::taxa_are_rows(physeq)) {
        otu_mat <- t(otu_mat)
    }
    
    # --- filter samples by depth ---
    keep <- rowSums(otu_mat) >= depth
    
    if (length(keep) != nrow(otu_mat)) {
        stop("Logical subsetting mismatch: check OTU table orientation.")
    }
    
    otu_mat <- otu_mat[keep, , drop = FALSE]
    
    n_samp <- nrow(otu_mat)
    n_taxa <- ncol(otu_mat)
    
    # --- single iteration ---
    if (n_iter == 1) {
        rare <- vegan::rrarefy(otu_mat, sample = depth)
        return(as.data.frame(rare))
    }
    
    # --- multiple iterations ---
    if (.array) {
        arr <- array(
            0,
            dim = c(n_samp, n_taxa, n_iter),
            dimnames = list(
                rownames(otu_mat),
                colnames(otu_mat),
                paste0("iter_", seq_len(n_iter))
            )
        )
        
        for (i in seq_len(n_iter)) {
            arr[, , i] <- vegan::rrarefy(otu_mat, sample = depth)
        }
        
        return(arr)
        
    } else {
        out <- vector("list", n_iter)
        names(out) <- paste0("iter_", seq_len(n_iter))
        
        for (i in seq_len(n_iter)) {
            out[[i]] <- as.data.frame(
                vegan::rrarefy(otu_mat, sample = depth)
            )
        }
        
        return(out)
    }
}


