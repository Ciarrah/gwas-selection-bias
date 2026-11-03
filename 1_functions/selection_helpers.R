# =============================================================================
# selection_helper.R
# Shared functions for generating selection into the sample.
# Contains functions only; no code runs when the file is sourced.
# =============================================================================

#' Calibrate a probit intercept to a target selection rate
#'
#' Finds a0 such that mean(pnorm(a0 + lin)) equals the target rate.
#'
#' @param lin Linear predictor of selection (without intercept).
#' @param rate Target proportion selected (between 0 and 1).
#' @return Intercept a0.
calibrate_intercept <- function(lin, rate) {
  stopifnot(rate > 0, rate < 1)
  stats::uniroot(function(a) mean(stats::pnorm(a + lin)) - rate,
                 interval = c(-20, 20))$root
}

#' Draw a selection indicator from a probit model
#'
#' @param lin Linear predictor of selection (without intercept).
#' @param rate Target proportion selected.
#' @return Integer 0/1 vector with attribute "intercept".
draw_selection <- function(lin, rate) {
  a0 <- calibrate_intercept(lin, rate)
  S  <- stats::rbinom(length(lin), 1, stats::pnorm(a0 + lin))
  attr(S, "intercept") <- a0
  S
}
