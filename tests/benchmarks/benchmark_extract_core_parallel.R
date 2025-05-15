# Comprehensive benchmark comparing original vs parallel versions
# Date: 2025-05-14
# User: jibarozzo

library(phyloseq)
library(BRCore)
library(microbenchmark)
library(parallel)
library(tidyverse)
library(vegan)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
source(here::here("tests/benchmarks/extract_core_parallel.R"))

# Load the GlobalPatterns dataset
data(GlobalPatterns, package = "phyloseq")

# Add grouping variable
sample_data(GlobalPatterns)$TestGroup <- sample(
  c("A", "B", "C"),
  nsamples(GlobalPatterns),
  replace = TRUE
)

# Create datasets of different sizes
datasets <- list(
  small = prune_taxa(
    names(sort(taxa_sums(GlobalPatterns), decreasing = TRUE)[1:100]),
    GlobalPatterns
  ),
  medium = prune_taxa(
    names(sort(taxa_sums(GlobalPatterns), decreasing = TRUE)[1:500]),
    GlobalPatterns
  ),
  large = prune_taxa(
    names(sort(taxa_sums(GlobalPatterns), decreasing = TRUE)[1:5000]),
    GlobalPatterns
  )
)

# Define parameters for testing
original_params <- list(
  Var = "SampleType",
  method = "increase",
  increase_value = 2
)

parallel_params <- list(
  Var = "SampleType",
  method = "increase",
  increase_value = 2,
  .parallel = TRUE,
  ncores = 4
)

# Function to benchmark with a specific dataset
run_size_benchmark <- function(dataset_name) {
  dataset <- datasets[[dataset_name]]

  cat("Benchmarking", dataset_name, "dataset (", ntaxa(dataset), "OTUs)\n")

  # Set parameters
  params <- c(list(physeq = dataset), original_params)
  params_parallel <- c(list(physeq = dataset), parallel_params)

  # Run benchmark
  result <- tryCatch(
    {
      microbenchmark(
        original = do.call(extract_core, params),
        parallel = do.call(extract_core_parallel, params_parallel),
        times = 50
      )
    },
    error = function(e) {
      cat("Error with", dataset_name, "dataset:", e$message, "\n")
      return(NULL)
    }
  )

  return(result)
}

# Function to run all benchmarks
run_all_benchmarks <- function() {
  results <- list()

  for (dataset_name in names(datasets)) {
    cat(
      "Starting benchmark for",
      dataset_name,
      "at",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "\n"
    )
    results[[dataset_name]] <- run_size_benchmark(dataset_name)
    cat(
      "Completed benchmark for",
      dataset_name,
      "at",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "\n"
    )
  }

  return(results)
}

# Main benchmark
cat("Starting all benchmarks at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
benchmark_results <- run_all_benchmarks()
cat(
  "All benchmarks completed at",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "\n"
)

# Process and save results
results_summary <- data.frame()

for (dataset_name in names(benchmark_results)) {
  if (!is.null(benchmark_results[[dataset_name]])) {
    summary_data <- summary(benchmark_results[[dataset_name]])
    summary_data$sd <- sd(benchmark_results[[dataset_name]]$time) / 1e9 # Convert to nanoseconds
    summary_data$dataset <- dataset_name
    summary_data$otu_count <- ntaxa(datasets[[dataset_name]])
    results_summary <- rbind(results_summary, summary_data)

    # Append to results summary
    results_summary <- rbind(results_summary, summary_data)
  }
}

# Save results
write.csv(
  results_summary,
  here::here("tests/benchmarks/extract_core_benchmark_summary.csv"),
  row.names = FALSE
)
saveRDS(
  benchmark_results,
  here::here("extract_core_benchmark_results_by_size.rds")
)

# Create summary plot
if (nrow(results_summary) > 0) {
  library(ggplot2)

  p <- ggplot(results_summary, aes(x = dataset, y = mean, fill = expr)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(
      aes(ymin = mean - sd, ymax = mean + sd),
      position = position_dodge(0.9),
      width = 0.2
    ) +
    theme_minimal() +
    labs(
      title = "extract_core Benchmark by Dataset Size",
      subtitle = paste("Benchmark run on", format(Sys.Date(), "%Y-%m-%d")),
      x = "Dataset Size",
      y = "Time (nanoseconds)",
      fill = "Implementation"
    ) +
    scale_y_log10()

  p
  # Save plot
  ggsave(
    here::here("tests/benchmarks/extract_core_benchmark_by_size.png"),
    p,
    width = 10,
    height = 6
  )
}

cat("Results processed and saved\n")
