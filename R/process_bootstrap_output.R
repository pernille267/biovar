#' @title Process and Summarize Bootstrap Model Output
#'
#' @param fit An \code{list}. The output from \code{bv_anova_bootstrap_ci()}
#' @param log_transformed A \code{logical} value. If \code{TRUE}, it assumes
#'   the model was fitted on log-transformed data, and results
#'   (like \code{beta}) are back-transformed. CV calculations are also adjusted
#'   accordingly.
#' @param analyte A \code{character} string used to label the measurand in the
#'   output tables. Defaults to "X".
#' @param group A \code{character} string used to label the group or cohort in
#'   the output tables. Defaults to "Y".
#' @param data An optional \code{data.frame} or \code{data.table} containing the
#'   original input data. If provided, its \code{SubjectID} column will be used
#'   to label the subjects in the subject-specific output table.
#'   Defaults to \code{NULL}.
#' @param digits An \code{integer} specifying the number of decimal places to use
#'   when formatting the confidence intervals. Defaults to 1.
#'
#' @returns A \code{data.table}.
#' @export
#'
#' @examples print(1)
process_bootstrap_output <- function(fit, log_transformed, analyte = "X", group = "Y", data = NULL, digits = 1L){
  
  # Format CI
  format_ci <- function(x, m) {
    if (!is.numeric(x) || length(x) != 3) {
      stop("Input 'x' must be a numeric vector of length 3.")
    }
    if (!is.numeric(m) || length(m) != 1 || m < 0 || m %% 1 != 0) {
      stop("Input 'm' must be a single non-negative integer.")
    }
    count_decimals <- function(val) {
      if ((val %% 1) != 0) {
        nchar(strsplit(as.character(val), ".", fixed = TRUE)[[1]][2])
      } else {
        0
      }
    }
    max_actual_decimals <- max(sapply(x, count_decimals), na.rm = TRUE)
    digits_to_show <- min(max_actual_decimals, m)
    format_string <- paste0("%.", digits_to_show, "f")
    result <- sprintf(
      "%s (%s-%s)",
      sprintf(format_string, x[1]),
      sprintf(format_string, x[2]),
      sprintf(format_string, x[3])
    )
    return(result)
  }
  
  beta <- c(fit$point_estimates$beta, fit$conf_intervals$beta_CI)
  SD_I <- c(fit$point_estimates$sigma_I, fit$conf_intervals$sigma_I_CI)
  SD_A <- c(fit$point_estimates$sigma_A, fit$conf_intervals$sigma_A_CI)
  SD_G <- c(fit$point_estimates$sigma_G, fit$conf_intervals$sigma_G_CI)
  HBHR <- c(fit$point_estimates$HBHR, fit$conf_intervals$HBHR_CI)
  
  CV_I <- NA
  CV_A <- NA
  CV_G <- NA
  
  if (log_transformed) {
    beta <- exp(beta)
    CV_I <- sqrt(exp(SD_I^2) - 1) * 100
    CV_A <- sqrt(exp(SD_A^2) - 1) * 100
    CV_G <- sqrt(exp(SD_G^2) - 1) * 100
  }
  
  else {
    CV_I <- SD_I / beta[1] * 100 
    CV_A <- SD_A / beta[1] * 100
    CV_G <- SD_G / beta[1] * 100
  }
  
  n_subjects <- NA_integer_
  n_results <- NA_integer_
  
  if (!is.null(data)) {
    n_subjects <- length(unique(data$SubjectID))
    n_results <- length(unique(paste0(data$SubjectID,"-",data$SampleID)))
  }
  
  output_1 <- data.table(
    "Measurand" = analyte,
    "Group" = group,
    "No. of individuals" = n_subjects,
    "No. of results" = n_results,
    "Mean concentration (95% CI)" = format_ci(beta, m = digits),
    "CVI (95% CI) %" = format_ci(CV_I, m = digits),
    "dCV_P(i) 50% (20%-80%) %" = NA,
    "CVG (95% CI) %" = format_ci(CV_G, m = digits),
    "CVA (95% CI) %" = format_ci(CV_A, m = digits),
    "HBHR (95% CI)" = format_ci(HBHR, m = digits)
  )
  
  return(output_1)
  
}