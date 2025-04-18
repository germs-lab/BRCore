#' Extract Core Microbial Taxa
#'
#' Identifies core microbial taxa based on abundance-occupancy distributions and their contributions
#' to Bray-Curtis dissimilarity. Core taxa are selected using either the "increase" or "elbow" method.
#'
#' The functions below use abundance-occupancy distributions fitted to a neutral model
#' described by Shade and Stopnisek, 2019. Core microbial taxa are selected based on their
#' contributions to overall microbial beta-diversity. In the described function, a core microbial taxa
#' must contribute at least ~2% of variation to Bray-Curtis dissimilarity to be considered a 'core' microbial taxa but this value can be manipulated within the function.

#' Core microbial taxa whose abundance and occupancy are above the fitted neutral model's confidence intervals
#' indicates greater occupancy across samples given abundance, suggesting deterministic selection by the plant.
#' Alternatively, core microbial taxa that are below the fitted neutral model indicates greater abundance given
#' lower occupancy; these taxa may be  dispersal limited.

#' More info on the functions can be found in Shade and Stopnisek, 2019. Original functions were developed in VanWallendael et al 2021.
#' Code was adapted and updated by Nicco Benucci (GLBRC, Bonito Lab) and Brandon Kristy (GLBRC, Evans Lab).
#' Error handling, helper functions and refactoring by Bolívar Aponte Rolón (CABBI, GERMS Lab) and Brandon Kristy. -2025-03-06
#' NOTE: This technique requires en even sampling depth to perform, which requires rarefaction.

#'
#' @param physeq A `phyloseq` object containing an OTU table, taxonomy table, and sample metadata.
#' @param Var A character string specifying the column name in the sample metadata to group samples.
#' @param method A character string specifying the method for selecting core taxa.
#' @param increase_value A numeric value specifying the threshold for core taxa selection.
#' @param Group A character string specifying the column name in the sample metadata to subset the data.
#' @param Level A character vector specifying the level(s) within the `Group` column to subset the data.
#'
#' @return A list containing the following elements:
#'   - `core_otus`: A vector of core OTU IDs.
#'   - `bray_curtis_ranked`: A data frame of Bray-Curtis dissimilarity rankings.
#'   - `otu_rankings`: A data frame of OTU rankings based on occupancy and abundance.
#'   - `occupancy_abundance`: A data frame of occupancy and abundance values.
#'   - `otu_table`: The OTU table used for analysis.
#'   - `sample_metadata`: The sample metadata used for analysis.
#'   - `taxonomy_table`: The taxonomy table used for analysis.
#'
#' @import phyloseq
#' @import dplyr
#' @import tidyr
#' @import vegan
#' @import tibble
#'
#' @examples
#' # Load the esophagus dataset
#' library(phyloseq)
#' data(esophagus, package = "phyloseq")
#'
#' #' Get the taxa names from the esophagus dataset
#' taxa_names <- taxa_names(esophagus)
#'
#' # Define realistic taxonomy levels
#' kingdoms <- c("Bacteria", "Archaea")
#' phyla <- c("Firmicutes",
#' "Bacteroidetes",
#' "Proteobacteria",
#' "Actinobacteria",
#' "Euryarchaeota")
#' classes <- c("Clostridia",
#' "Bacteroidia",
#' "Gammaproteobacteria",
#' "Actinobacteria",
#' "Methanobacteria")
#' orders <- c("Clostridiales",
#' "Bacteroidales",
#' "Enterobacterales",
#' "Bifidobacteriales",
#' "Methanobacteriales")
#' families <- c("Lachnospiraceae",
#' "Bacteroidaceae",
#' "Enterobacteriaceae",
#' "Bifidobacteriaceae",
#' "Methanobacteriaceae")
#' genera <- c("Blautia",
#' "Bacteroides",
#' "Escherichia",
#' "Bifidobacterium",
#' "Methanobrevibacter")
#' species <- c("Blautia producta",
#' "Bacteroides fragilis",
#' "Escherichia coli",
#' "Bifidobacterium longum",
#' "Methanobrevibacter smithii")
#'
#' # Create a mock taxonomy table
#' set.seed(8998); mock_taxonomy_table <- matrix(
#'   c(
#'     sample(kingdoms, length(taxa_names), replace = TRUE),
#'     sample(phyla, length(taxa_names), replace = TRUE),
#'     sample(classes, length(taxa_names), replace = TRUE),
#'     sample(orders, length(taxa_names), replace = TRUE),
#'     sample(families, length(taxa_names), replace = TRUE),
#'     sample(genera, length(taxa_names), replace = TRUE),
#'     sample(species, length(taxa_names), replace = TRUE)
#'   ),
#'   nrow = length(taxa_names),
#'   ncol = 7,
#'   dimnames = list(
#'     taxa_names,
#'     c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#'   )
#' )
#'
#' # Convert to a taxonomy table object
#' tax_table <- tax_table(mock_taxonomy_table)
#'
#' # Add sample metadata to the esophagus dataset
#' sample_data <- data.frame(
#'   Sample = sample_names(esophagus),
#'   Group = sample(c("A", "B"), nsamples(esophagus), replace = TRUE), #' Random groups
#'   row.names = sample_names(esophagus)
#' )
#'
#' # Add the taxonomy table to the esophagus dataset
#' esophagus_with_tax <- merge_phyloseq(esophagus, tax_table, sample_data(sample_data))
#'
#' # Extract core taxa using the "increase" method
#' core_result <- extract_core(
#'   physeq = esophagus_with_tax,
#'   Var = "Group",
#'   method = "increase",
#'   increase_value = 2
#' )
#'
#' # View the results
#' print(core_result)
#' @export


extract_core <- function(physeq,
                         Var,
                         method,
                         increase_value = NULL,
                         Group = NULL,
                         Level = NULL) {
  set.seed(37920)
  
  # Error handling: type check
  if (!inherits(physeq, "phyloseq")) {
    cli::cli_abort(
      "{.arg physeq} must be a 'phyloseq' object.\nYou've supplied a {class(physeq)[1]} vector."
    )
  }
  # If the check passes, continue processing
  cli::cli_alert_success("Input phyloseq object is valid!")
  
  
  #################
  ## Rarefaction ##
  #################
  
  # input dataset needs to be rarified and minimum depth included
  nReads <- min(sample_sums(physeq))
  
  if (min(sample_sums(physeq)) == max(sample_sums(physeq))) {
    # nReads %T>% print()
    # rarefied %T>% print()
    taxon <- tax_table(physeq) %>%
      as.data.frame.matrix()
    
    dim(taxon) %T>% print()
  } else {
    nReads <- min(sample_sums(physeq))
    # nReads %T>% print()
    rarefied <-
      rarefy_even_depth(
        physeq,
        sample.size = nReads,
        trimOTUs = TRUE,
        replace = TRUE,
        verbose = FALSE
      )
    
    taxon <- tax_table(rarefied) %>%
      as.data.frame.matrix()
  }
  
  
  #################
  ## Subsetting ##
  #################
  
  # choosing a subset or using the whole phyloseq object as is
  if (is.null(Group)) {
    otu <- rarefied@otu_table %>% as("matrix")
    map <- rarefied@sam_data %>% as("data.frame")
  } else {
    sub_group <- substitute(Group)
    sub_set <- subset(sample_data(rarefied), eval(parse(text = sub_group)) %in% Level)
    physeq1 <- merge_phyloseq(otu_table(rarefied),
                              tax_table(rarefied),
                              refseq(rarefied),
                              sub_set)
    
    otu_table(physeq1) <- otu_table(physeq1)[which(rowSums(otu_table(physeq1)) > 0), ]
    otu <- physeq1@otu_table %>% as("matrix")
    map <- physeq1@sam_data %>% as("data.frame")
    print("Grouping Factor")
    map[, Group] %T>% print()
    
    taxon <- tax_table(physeq1) %>%
      as.data.frame.matrix()
  }
  
  map$SampleID <- rownames(map)

  
  #############################
  ## occupancy and abundance ##
  #############################
  
  # calculating occupancy and abundance
  otu_PA <-
    1 * ((otu > 0) == 1) # presence-absence data
  otu_occ <-
    rowSums(otu_PA) / ncol(otu_PA) # occupancy calculation
  otu_rel <-
    apply(decostand(otu, method = "total", MARGIN = 2), 1, mean) # mean relative abundance
  # tibble::rownames_to_column appears to create an error at line #106 for left_join(), sticking with add_rownames for now
  #   occ_abun <-
  #     add_rownames(as.data.frame(cbind(otu_occ, otu_rel)), "otu") # combining occupancy and abundance data frame
  # NOTE! add_rownames is deprecated and generates a warning, a bug of tidyverse,
  # alternative you can use:
  occ_abun <- tibble::rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)), "otu")
  
  
  ##################
  ## Ranking OTUs ##
  ##################
  
  # Ranking OTUs based on their occupancy
  # For calculating ranking index we included following conditions:
  # - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
  # - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)
  
  Var <- rlang::enquo(Var) # lazy evaluation
  PresenceSum <-
    data.frame(otu = as.factor(row.names(otu)), otu) %>%
    gather(SampleID, abun, -otu) %>%
    left_join(map, by = "SampleID") %>%
    group_by(otu, !!Var) %>%
    dplyr::summarise(
      time_freq = sum(abun > 0) / length(abun),
      # frequency of detection between time points
      coreTime = ifelse(time_freq == 1, 1, 0)
    ) %>% # 1 only if occupancy 1 with specific time, 0 if not
    group_by(otu) %>%
    dplyr::summarise(
      sumF = sum(time_freq),
      sumG = sum(coreTime),
      nS = length(!!Var) * 2,
      Index = (sumF + sumG) / nS
    ) # calculating weighting Index based on number of points detected
  
  
  # Ranked OTUs
  otu_ranked <- occ_abun %>%
    left_join(PresenceSum, by = "otu") %>%
    transmute(otu = otu, rank = Index) %>%
    arrange(desc(rank))

  
  ###############################
  ## Bray-Curtis Dissimilarity ##
  ###############################
  
  # Helper function: Calculate Bray-Curtis values
  calculate_bc <- function(matrix, nReads) {
    if (nrow(matrix) == 0) {
      cli::cli_alert_warning("{.arg matrix} is empty. Enter a non-empty matrix.")
      return(list(values = numeric(0), names = character(0)))
    }
    
    bc_values <- apply(combn(ncol(matrix), 2), 2, function(x) {
      sum(abs(matrix[, x[1]] - matrix[, x[2]])) / (2 * nReads)
    })
    
    x_names <- apply(combn(ncol(matrix), 2), 2, function(x) {
      paste(colnames(matrix)[x], collapse = "-")
    })
    
    list(values = bc_values, names = x_names)
  }
  
  # Calculating BC dissimilarity based on the 1st ranked OTU
  cli::cli_alert_info("Calculating BC dissimilarity based on the 1st ranked OTU")
  
  start_matrix <- t(as.matrix(otu[otu_ranked$otu[1], ]))
  first_bc <- calculate_bc(start_matrix, nReads)
  BCaddition <- data.frame(x_names = first_bc$names, "1" = first_bc$values)
  
  cli::cli_alert_success("BC dissimilarity based on the 1st ranked OTU complete")
  
  # Calculating BC dissimilarity based on additon of ranked OTUs from 2nd to nth.
  # Set to the entire length of OTUs in the dataset. It might take some time if more than 5000 OTUs are included.
  
  cli::cli_alert_info("Calculating BC dissimilarity based on ranked OTUs, starting at {Sys.time()}")
  progressbar_calc_bc <- cli::cli_progress_bar(
    name = "Calculating BC rankings",
    total = nrow(otu_ranked) - 1,
    format = "{cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
    .auto_close = TRUE,
    .envir = parent.frame()
  )
  
  for (i in 2:nrow(otu_ranked)) {
    # Add next OTU
    current_matrix <- rbind(start_matrix, t(otu[otu_ranked$otu[i], ]))
    
    # Calculate BC
    current_bc <- calculate_bc(current_matrix, nReads)
    
    # Update data frame
    df_a <- data.frame(x_names = current_bc$names, val = current_bc$values)
    names(df_a)[2] <- i
    BCaddition <- left_join(BCaddition, df_a, by = "x_names")
    
    cli::cli_progress_update(id = progressbar_calc_bc)
  }
  
  cli::cli_progress_done(id = progressbar_calc_bc)
  
  cli::cli_alert_success("BC ranks done!")
  
  
  
  # Ranking the obtained BC
  rownames(BCaddition) <- BCaddition$x_names
  temp_BC <- BCaddition
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)
  BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>%
    gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    summarise(MeanBC = mean(BC)) %>% # Calculate mean Bray-Curtis dissimilarity
    arrange(desc(-MeanBC)) %>%
    mutate(proportionBC = MeanBC / max(MeanBC)) # Calculate proportion of the dissimilarity explained by the n number of ranked OTUs

  
  # Calculating the increase BC
  Increase <- BC_ranked$MeanBC[-1] / BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC = c(0, (Increase)), rank = factor(c(1:(
    length(Increase) + 1
  ))))
  BC_ranked <- left_join(BC_ranked, increaseDF, by = join_by(rank))
  BC_ranked <- BC_ranked[-nrow(BC_ranked), ]
  BC_ranked <- drop_na(BC_ranked)
  
  
  #############################
  ## Increase in Bray-Curtis ##
  #############################
  
  # Calculating the increase BC
  Increase <- BC_ranked$MeanBC[-1] / BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC = c(0, (Increase)), rank = factor(c(1:(
    length(Increase) + 1
  ))))
  BC_ranked <- left_join(BC_ranked, increaseDF, by = join_by(rank, IncreaseBC))
  BC_ranked <- BC_ranked[-nrow(BC_ranked), ]
  BC_ranked <- drop_na(BC_ranked)
  
  
  # Error handling: Make sure a method is specified as 'increase' or 'elbow'
  # Check if method is provided and is one of the allowed values
  if (missing(method) || !method %in% c("increase", "elbow")) {
    cli::cli_abort("{.arg method} must be specified as either 'increase' or 'elbow'.",
                   call. = FALSE)
  }
  
  if (method == "elbow") {
    cli::cli_alert_success("Performing method 'elbow'")
    fo_difference <- function(pos) {
      left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
      right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
      core_otus <- otu_ranked$otu[1:elbow]
      # core_otus %T>% print()
    }
    occ_abun$fill <- "no"
    occ_abun$fill[occ_abun$otu %in% core_otus] <- "core"
  }
  
  # Creating threshold for core inclusion - last call method using a
  # final increase in BC similarity of equal or greater than increase_value
  if (method == "increase") {
    # Error handling: If method 'increase' is chosen, make sure that                                      'increase_value' is specified
    if (missing(increase_value) || is.null(increase_value)) {
      cli::cli_abort("{.arg increase_value} must be specified when method is 'increase'.")
    }
    if (!is.numeric(increase_value)) {
      cli::cli_abort("{arg increase_value} must be a numeric value.")
    }
    
    cli::cli_alert_success("Performing method 'increase'")
    
    # Convert the % increase into a decimal value and add 1
    perc_increase <- 1 + (increase_value * 0.01)
    
    lastCall <-
      as.numeric(as.character(
        dplyr::filter(BC_ranked, IncreaseBC >= perc_increase)$rank
      ))
    
    # Convert lastCall to character to subset core_otus by rowname
    lastCall <- as.character(lastCall)
    core_otus <- otu_ranked[rownames(otu_ranked) %in% lastCall, ]
    
    core_otus <- as.vector(core_otus$otu)
    
    # Filter non-core and core OTUs
    occ_abun$fill <- "no"
    occ_abun$fill[occ_abun$otu %in% core_otus] <- "core"
  }
  
  # Create named return list
  return_list <- list(
    core_otus = core_otus,
    bray_curtis_ranked = BC_ranked,
    otu_rankings = otu_ranked,
    occupancy_abundance = occ_abun,
    otu_table = otu,
    sample_metadata = map,
    taxonomy_table = taxon
  )
  
  # Validate return list
  if (any(sapply(return_list, is.null))) {
    cli::cli_alert_warning("One or more return list elements are NULL.")
  }
  
  return(return_list)
}
