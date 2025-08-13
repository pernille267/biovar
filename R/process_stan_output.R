#' Process and Summarize Stan Model Output
#'
#'
#' @param fit An object, typically a \code{list} generated from
#'   \code{rstan::extract}, containing the posterior draws from the Stan model.
#'   It is expected to contain vectors and matrices of draws for parameters
#'   like \code{beta}, \code{sigma_I}, \code{sigma_A}, \code{sigma_G},
#'   \code{df_I}, \code{df_A}, etc.
#' @param log_transformed A \code{logical} value. If \code{TRUE}, it assumes the model
#'   was fitted on log-transformed data, and results (like \code{beta}) are
#'   back-transformed. CV calculations are also adjusted accordingly.
#' @param analyte A \code{character} string used to label the measurand in the
#'   output tables. Defaults to "X".
#' @param group A \code{character} string used to label the group or cohort in
#'   the output tables. Defaults to "Y".
#' @param data An optional \code{data.frame} or \code{data.table} containing the
#'   original input data. If provided, its \code{SubjectID} column will be used to
#'   label the subjects in the subject-specific output table.
#'   Defaults to \code{NULL}.
#' @param digits An \code{integer} specifying the number of decimal places to use
#'   when formatting the credible intervals. Defaults to 1.
#' 
#' @description
#' Takes the posterior draws from a Stan model fit and computes key biological
#' and analytical metrics, such as coefficients of variation (CVs) and credible
#' intervals, formatting them into user-friendly tables.
#'
#' @details
#' This function serves as the primary post-processing step. It calculates the
#' mean concentration, within-subject CV (\code{CVI}), between-subject CV (\code{CVG}),
#' and analytical CV (\code{CVA}) from the posterior distributions.
#' It handles calculations differently depending on whether the model was run on a
#' log-transformed scale. The function also computes the distribution of the
#' predictive within-subject CV (\code{dCV_P(i)}) and the heterogeneity of this
#' distribution (HBHR).
#'
#' \strong{Important Note on Model Validity:}
#' The reliability of the metrics produced by this function depends entirely
#' on the quality and convergence of the Stan model fit. The \code{fit} object is
#' typically the result of \code{rstan::extract(your_stanfit_object)}. Before
#' trusting the output of this function, it is essential to perform standard MCMC
#' diagnostics. Key checks include:
#' \itemize{
#'   \item Verifying that there are no divergent transitions after warmup.
#'   \item Checking for transitions that exceeded the maximum treedepth.
#'   \item Ensuring that the R-hat (potential scale reduction factor) for all parameters is close to 1.0 (e.g., < 1.01).
#'   \item Confirming that the effective sample size (\code{n_eff}) is sufficiently large for all parameters of interest.
#' }
#' These diagnostics can be checked by inspecting the summary of the original \code{stanfit} object.
#'
#' @return A \code{list} containing three \code{data.table} objects:
#' \describe{
#'   \item{Table 1}{A single-row summary table containing the overall model
#'     estimates, including mean concentration, CVI, CVG, CVA, the distribution
#'     of the predictive within-subject CV, and the HBHR, all with 95\% credible
#'     intervals.}
#'   \item{Table 2}{A subject-specific table showing the individual within-subject
#'     CV (\code{CV_P(i)}) for each subject, including the mean, median, and 95\%
#'     credible interval. Subjects are ordered by their median \code{CV_P(i)}.}
#'   \item{Table 3}{A summary table for the model's degrees of freedom
#'     parameters (\code{df_I} and \code{df_A}), showing their mean values and 95\%
#'     credible intervals.}
#' }
#'
#' @seealso \code{\link{plot_subject_specific_CVI}} for visualizing the results.
#' @export

process_stan_output <- function(fit, log_transformed, analyte = "X", group = "Y", data = NULL, digits = 1L) {
  
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
  
  beta <- NA_real_
  CV_I <- NA_real_
  CV_I_pred <- NA_real_
  CV_pi <- NA_real_
  CV_A <- NA_real_
  CV_G <- NA_real_
  HBHR <- NA_real_
  HBHR_pred <- NA_real_
  
  if (log_transformed) {
    beta <- exp(fit$beta)
    CV_I <- sqrt(exp(fit$sigma_I^2 * (mean(fit$df_I) - 2.0) / mean(fit$df_I)) - 1) * 100
    CV_A <- sqrt(exp(fit$sigma_A^2 * (mean(fit$df_A) - 2.0) / mean(fit$df_A)) - 1) * 100
    CV_G <- sqrt(exp(fit$sigma_G^2) - 1) * 100
    CV_I_pred <- sqrt(exp(fit$sigma_I_pred^2) - 1) * 100
    CV_pi <- sqrt(exp(fit$sigma_i^2 * (mean(fit$df_I) - 2.0) / mean(fit$df_I)) - 1) * 100
  }
  
  else {
    beta <- fit$beta
    CV_I <- fit$sigma_I * sqrt((mean(fit$df_I) - 2) / mean(fit$df_I)) / mean(fit$beta) * 100 
    CV_A <- fit$sigma_A * sqrt((mean(fit$df_A) - 2) / mean(fit$df_A)) / mean(fit$beta) * 100
    CV_G <- fit$sigma_G / mean(fit$beta) * 100
    CV_I_pred <- fit$sigma_I_pred / mean(fit$beta) * 100 
    CV_pi <- fit$sigma_i * sqrt((mean(fit$df_I) - 2) / mean(fit$df_I)) / mean(fit$beta) * 100
  }
  
  
  output_1 <- data.table(
    "Measurand" = analyte,
    "Group" = group,
    "No. of individuals" = ncol(fit$G_raw),
    "No. of results" = ncol(fit$I_raw),
    "Mean concentration (95% CrI)" = format_ci(c(mean(beta), quantile(beta, probs = 0.025), quantile(beta, probs = 0.975)), m = digits),
    "CVI (95% CrI) %" = format_ci(c(mean(CV_I, na.rm = TRUE), quantile(CV_I, probs = 0.025, na.rm = TRUE), quantile(CV_I, probs = 0.975, na.rm = TRUE)), m = digits),
    "dCV_P(i) 50% (20%-80%) %" = format_ci(quantile(CV_I_pred, probs = c(0.50, 0.20, 0.80)), m = digits),
    "CVG (95% CrI) %" = format_ci(c(mean(CV_G, na.rm = TRUE), quantile(CV_G, probs = 0.025, na.rm = TRUE), quantile(CV_G, probs = 0.975, na.rm = TRUE)), m = digits),
    "CVA (95% CrI) %" = format_ci(c(mean(CV_A, na.rm = TRUE), quantile(CV_A, probs = 0.025, na.rm = TRUE), quantile(CV_A, probs = 0.975, na.rm = TRUE)), m = digits),
    "HBHR %" = format(sd(CV_I_pred) / mean(CV_I_pred) * 100, nsmall = 1, digits = 1)
  )
  
  subjects <- data.table("SubjectID" = 1:ncol(CV_pi))
  output_2 <- apply(
    X = CV_pi,
    MARGIN = 2,
    FUN = function(SUBJECT) {
      list("mean_CV_P(i)" = mean(SUBJECT, na.rm = TRUE),
           "median_CV_P(i)" = median(SUBJECT, na.rm = TRUE),
           "CV_P(i)_lwr" = quantile(SUBJECT, probs = 0.025, na.rm = TRUE),
           "CV_P(i)_upr" = quantile(SUBJECT, probs = 0.975, na.rm = TRUE),
           "dCV_P(i)_median" = median(CV_I_pred, na.rm = TRUE),
           "dCV_P(i)_lwr" = quantile(CV_I_pred, probs = 0.20, na.rm = TRUE),
           "dCV_P(i)_upr" = quantile(CV_I_pred, probs = 0.80, na.rm = TRUE))
    }
  )
  
  output_2 <- rbindlist(output_2)
  output_2 <- cbind(subjects, output_2)
  
  if (!is.null(data)) {
    output_2$SubjectID <- unique(data$SubjectID)
  }
  
  output_2$SubjectID <- reorder(x = output_2$SubjectID, X = output_2$`median_CV_P(i)`)
  
  output_3 <- data.table(
    "Measurand" = analyte,
    "Group" = group,
    "No. of individuals" = ncol(fit$G_raw),
    "No. of results" = ncol(fit$I_raw),
    "Mean DFI (95% CrI)" = format_ci(c(mean(fit$df_I, na.rm = TRUE), quantile(fit$df_I, probs = 0.025, na.rm = TRUE), quantile(fit$df_I, probs = 0.975, na.rm = TRUE)), m = digits),
    "Mean DFA (95% CrI)" = format_ci(c(mean(fit$df_A, na.rm = TRUE), quantile(fit$df_A, probs = 0.025, na.rm = TRUE), quantile(fit$df_A, probs = 0.975, na.rm = TRUE)), m = digits)
  )
  
  return(list(output_1, output_2, output_3))
  
}

#' Plot Subject-Specific Within-Subject Coefficients of Variation
#'
#' @description
#' Generates a caterpillar plot visualizing the median and 95% credible intervals
#' for the subject-specific within-subject CV (\code{CVp(i)}).
#'
#' @details
#' This function first calls \code{\link{process_stan_output}} to get the
#' processed data. It then uses the subject-specific results (the second element
#' of the returned list) to create a \code{ggplot2} object.
#'
#' The plot displays:
#' \itemize{
#'   \item A point for each subject's median \code{CVp(i)}.
#'   \item Horizontal error bars representing the 95% credible interval for each subject's \code{CVp(i)}.
#'   \item Vertical red lines indicating the population-level distribution of the
#'     predictive within-subject CV (solid line for the median, dashed lines for the
#'     20th and 80th percentiles). In other words, the 60% prediction interval.
#' }
#' Subjects are ordered on the y-axis by their median \code{CVp(i)}.
#'
#' @param fit An object containing the posterior draws from the Stan model. Passed directly to \code{\link{process_stan_output}}.
#' @param log_transformed A \code{logical} value indicating if the model was fitted on log-transformed data. Passed directly to \code{\link{process_stan_output}}.
#' @param analyte A \code{character} string for the measurand label. Passed directly to \code{\link{process_stan_output}}.
#' @param group A \code{character} string for the group label. Passed directly to \code{\link{process_stan_output}}.
#' @param data An optional \code{data.frame} or \code{data.table} with original subject IDs. Passed directly to \code{\link{process_stan_output}}.
#' @param digits An \code{integer} for rounding. Passed directly to \code{\link{process_stan_output}}.
#'
#' @return A \code{ggplot2} object representing the plot. This can be further
#'   customized by the user (e.g., by adding themes or other layers).
#'
#' @seealso \code{\link{process_stan_output}} which generates the data for this plot.
#' @import ggplot2
#' @importFrom stats reorder
#' @export
#' @examples
#' # This is a mock example as it requires a real Stan fit object.
#' # 1. Create a fake 'fit' object with the required structure
#' n_subjects <- 10
#' n_draws <- 1000
#' fake_fit <- list(
#'   beta = rnorm(n_draws, 100, 5),
#'   sigma_I = rlnorm(n_draws, log(10), 0.1),
#'   sigma_A = rlnorm(n_draws, log(3), 0.1),
#'   sigma_G = rlnorm(n_draws, log(15), 0.1),
#'   sigma_I_pred = rlnorm(n_draws, log(10), 0.2),
#'   sigma_i = matrix(rlnorm(n_draws * n_subjects, log(10), 0.3), nrow = n_draws),
#'   df_I = rep(30, n_draws),
#'   df_A = rep(50, n_draws),
#'   G_raw = matrix(0, nrow = 1, ncol = n_subjects), # For counting
#'   I_raw = matrix(0, nrow = 1, ncol = n_subjects * 3) # For counting
#' )
#'
#' # 2. Generate the plot
#' p <- plot_subject_specific_CVI(
#'   fit = fake_fit,
#'   log_transformed = FALSE,
#'   analyte = "Testosterone",
#'   group = "Healthy Adults"
#' )
#'
#' # print(p) # Uncomment to display the plot
#'

plot_subject_specific_CVI <- function(fit, log_transformed, analyte = "X", group = "Y", data = NULL, digits = 1L) {
  processed_stan_output <- process_stan_output(fit = fit,
                                               log_transformed = log_transformed,
                                               analyte = analyte,
                                               group = group,
                                               data = data,
                                               digits = digits)
  
  plotting_data <- processed_stan_output[[2]]
  
  ggplot(data = plotting_data) +
    geom_vline(xintercept = c(plotting_data$`dCV_P(i)_lwr`[1], plotting_data$`dCV_P(i)_upr`[1]),
               color = "red",
               linewidth = 0.75,
               linetype = 2) +
    geom_vline(xintercept = plotting_data$`dCV_P(i)_median`,
               color = "red",
               linewidth = 0.75) +
    geom_errorbarh(mapping = aes(y = SubjectID, xmin = `CV_P(i)_lwr`, xmax = `CV_P(i)_upr`),
                   linewidth = 0.75) +
    geom_point(mapping = aes(y = SubjectID, x = `median_CV_P(i)`),
               size = 3,
               shape = 21,
               fill = "orange") +
    scale_x_continuous(name = expression(CV[p(i)]),
                       n.breaks = 10,
                       labels = function(x) paste(format(x, nsmall = 1, digits = 1), " %")) +
    ylab(label = "Subject") +
    labs(title = "Subject-specific CVs with 95% CrI",
         subtitle = paste0("Measurand: ", analyte, ", ", "Cohort: ", group)) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "orange"),
          plot.subtitle = element_text(face = "bold", hjust = 0.5, color = "black"))
  
}






