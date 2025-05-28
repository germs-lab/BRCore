#' Detect HPC environment type
#'
#' Detects if code is running in an HPC environment and identifies which scheduler.
#'
#' @return Character string identifying the HPC environment type:
#'   "slurm", "pbs", "sge", "lsf", "condor", "generic_hpc", or "local"
#'
#' @examples
#' detect_hpc_type()
#' 
#' @export
detect_hpc_type <- function() {
  # Map of environment variables to HPC types
  hpc_type_vars <- list(
    slurm = c("SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_JOB_NAME"),
    pbs = c("PBS_JOBID", "PBS_NUM_PPN", "PBS_O_WORKDIR"),
    sge = c("SGE_TASK_ID", "NSLOTS", "SGE_ROOT"),
    lsf = c("LSF_JOBID", "LSB_DJOB_NUMPROC", "LSB_JOBID"),
    condor = c("CONDOR_JOB_ID", "_CONDOR_JOB_AD", "CONDOR_SLOT")
  )

  # Generic HPC variables (not specific to a scheduler)
  generic_vars <- c("OMP_NUM_THREADS", "JOB_ID", "HOSTNAME")

  # Check each HPC type
  for (type in names(hpc_type_vars)) {
    if (
      any(sapply(hpc_type_vars[[type]], function(var) {
        val <- Sys.getenv(var, "")
        !identical(val, "")
      }))
    ) {
      return(type)
    }
  }

  # Check generic HPC indicators
  if (
    any(sapply(generic_vars, function(var) {
      val <- Sys.getenv(var, "")
      !identical(val, "")
    }))
  ) {
    return("generic_hpc")
  }

  # Not in an HPC environment
  return("local")
}

#' Check if running in an HPC environment
#'
#' @return Logical indicating whether we're in any HPC environment
#' @examples
#' is_hpc_environment()
#' 
#' @export
is_hpc_environment <- function() {
  !identical(detect_hpc_type(), "local")
}

#' Get the number of available CPU cores based on environment
#'
#' Determines the appropriate number of cores to use based on whether
#' we're in an HPC environment or a local machine.
#'
#' @param default Default number of cores to use if detection fails
#' @return Integer number of cores available
#' @examples
#' get_available_cores()
#' get_available_cores(default = 2)
#' 
#' @import parallelly
#' @import parallel
#' 
#' @export
get_available_cores <- function(default = 1) {
  if (is_hpc_environment()) {
    # Check HPC-specific environment variables in a prioritized order
    hpc_core_vars <- c(
      "SLURM_CPUS_PER_TASK", # SLURM
      "PBS_NUM_PPN", # PBS/TORQUE
      "NSLOTS", # SGE
      "LSB_DJOB_NUMPROC", # LSF
      "OMP_NUM_THREADS" # OpenMP (generic)
    )

    # Try each variable in order
    for (var in hpc_core_vars) {
      core_str <- Sys.getenv(var, "")
      if (core_str != "") {
        cores <- suppressWarnings(as.integer(core_str) - 1L) # Keeping 1 core available
        if (!is.na(cores) && cores > 0) {
          cli::cli_alert_info(
            "HPC environment detected. Using {cores} core(s)/worker(s) from {var}.
            Keeping 1 core/worker available for system."
          )
          return(cores)
        }
      }
    }

    # If no specific HPC variable found, try parallelly if available
    if (requireNamespace("parallelly", quietly = TRUE)) {
      cores <- parallelly::availableCores() - 1L
      cli::cli_alert_info(
        "HPC environment detected. Using {cores} core(s)/worker(s) via parallelly::availableCores().
        Keeping 1 core/worker available for system."
      )
      return(cores)
    }

    # Last resort for HPC: use a safe default
    cli::cli_alert_warning(
      "HPC environment detected but could not determine core count. Using {default} core(s)."
    )
    return(default)
  } else {
    # Local machine: use parallelly if available, otherwise safer detectCores
    if (requireNamespace("parallelly", quietly = TRUE)) {
      cores <- parallelly::availableCores() - 1L
      cli::cli_alert_info(
        "Local environment. Using {cores} core(s)/worker(s) via parallelly::availableCores().
        Keeping 1 core/worker available for system."
      )
      return(cores)
    } else {
      # Safe local default: leave one core free for system
      cores <- max(1, parallel::detectCores() - 1L)
      if (is.na(cores)) cores <- default
      cli::cli_alert_info(
        "Local environment. Using {cores} core(s)/worker(s) via parallel::detectCores().
        Keeping 1 core/worker available for system."
      )
      return(cores)
    }
  }
}

#' Set up an appropriate parallel backend
#'
#' Creates a parallel cluster if in an HPC environment or returns
#' the number of cores to use if on a local machine.
#'
#' @param default Default number of cores to use if detection fails
#' @return Either a cluster object (HPC) or integer number of cores (local)
#' @examples
#' setup_parallel_backend()
#' 
#' @import parallelly
#' @import parallel
#' @export
setup_parallel_backend <- function(default = 1) {
  if (is_hpc_environment()) {
    if (requireNamespace("parallelly", quietly = TRUE)) {
      # Create a PSOCK cluster with available workers
      workers <- parallelly::availableWorkers()
      cli::cli_alert_info(
        "Creating HPC cluster with {length(workers) - 1L} core(s)/worker(s) via parallely::makeClusterPSOCK().
        Keeping 1 core/worker available for system."
      )
      cl <- parallelly::makeClusterPSOCK(length(workers) - 1L, autoStop = TRUE) # Keep one for other processes
      return(cl)
    } else {
      # Fallback if parallelly is not available
      cores <- get_available_cores(default)
      cli::cli_alert_info(
        "Creating HPC cluster with {cores - 1L} core(s)/worker(s) using parallel::makeCluster(). 
        Keeping 1 core/worker available for system."
      )
      cl <- parallel::makeCluster(cores - 1L)
      return(cl)
    }
  } else {
    # Just return the number of cores for local processing
    cores <- get_available_cores(default)
    return(cores)
  }
}
