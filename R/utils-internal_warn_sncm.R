#' Internal warning handler for sncm.fit likelihood functions
#'
#' Wraps a function call with a `withCallingHandlers()` that captures warnings
#' and re-emits them via `cli::cli_warn()`. Deduplicates repeated warnings
#' across optimizer iterations using an external flag and setter.
#'
#' @param fun A zero-argument function (lambda) to evaluate.
#' @param message A string prefix for the warning message.
#' @param warned Logical flag indicating whether a warning has already fired.
#' @param warn_setter A zero-argument function that sets the external `warned`
#'   flag to `TRUE` via `<<-`.
#'
#' @return The result of `fun()`.
#' @keywords internal
.internal_warn_sncm <- function(
  fun,
  message = NULL,
  warned = FALSE,
  warn_setter = NULL
) {
  withCallingHandlers(
    fun(),
    warning = function(w) {
      if (!warned) {
        cli::cli_warn("{message}: {conditionMessage(w)}")
        if (!is.null(warn_setter)) warn_setter()
      }
      invokeRestart("muffleWarning")
    }
  )
}
