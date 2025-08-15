#' Format a credible interval from a vector of three values.
#' @param stats_vec A numeric vector of length 3: c(point_estimate, lower, upper).
#' @param digits The number of decimal places to format to.
#' @return A formatted string like "point_est (lower-upper)".
format_ci <- function(stats_vec, digits) {
  if (!is.numeric(stats_vec) || length(stats_vec) != 3) {
    # Return NA for invalid input to avoid crashing the whole summary
    return(NA_character_)
  }
  
  format_string <- paste0("%.", as.integer(digits), "f")
  
  sprintf(
    "%s (%s-%s)",
    sprintf(format_string, stats_vec[1]),
    sprintf(format_string, stats_vec[2]),
    sprintf(format_string, stats_vec[3])
  )
}

#' Calculate CVs assuming a log-transformed model.
#' @param fit A list of posterior draws.
#' @return A list containing posterior draws for each calculated CV.
calculate_cvs_log_scale <- function(fit) {
  # For robustness, use the mean of the df posteriors as a fixed value
  # This avoids a complex calculation for every single draw.
  mean_df_I <- mean(fit$df_I)
  mean_df_A <- mean(fit$df_A)
  
  # --- Define transformation functions for clarity ---
  
  # CV for a log-Student-t distribution
  logt_scale_to_cv <- function(scale, df) {
    # The variance of a log-t is equivalent to a log-normal with sd = sqrt(scale^2 * df / (df-2))
    sqrt(exp(scale^2 * df / (df - 2)) - 1) * 100
  }
  
  # CV for a log-normal distribution
  lognormal_sd_to_cv <- function(log_sd) {
    sqrt(exp(log_sd^2) - 1) * 100
  }
  
  # --- Apply transformations to posterior draws ---
  beta_pred <- exp(fit$beta)
  cv_I <- logt_scale_to_cv(fit$sigma_I, mean_df_I)
  cv_A <- logt_scale_to_cv(fit$sigma_A, mean_df_A)
  cv_G <- lognormal_sd_to_cv(fit$sigma_G)
  cv_I_pred <- lognormal_sd_to_cv(fit$sigma_I_pred)
  
  # For individual CVs (CV_pi), this is a matrix operation.
  # We apply the log-t transform to each column of fit$sigma_i.
  # Since fit$sigma_i is (draws x subjects), this is already vectorized.
  cv_pi <- logt_scale_to_cv(fit$sigma_i, mean_df_I)
  
  return(list(
    beta_pred = beta_pred, cv_I = cv_I, cv_A = cv_A, cv_G = cv_G,
    cv_I_pred = cv_I_pred, cv_pi = cv_pi
  ))
}

#' Calculate CVs assuming a model on the original (identity) scale.
#' @param fit A list of posterior draws.
#' @return A list containing posterior draws for each calculated CV.
calculate_cvs_identity_scale <- function(fit) {
  
  # Predictive Distribution of CV_i
  n_draws <- length(fit$beta)
  mu_new_subject <- rnorm(n = n_draws, mean = fit$beta, sd = fit$sigma_G)
  sd_new_subject <- fit$sigma_I_pred
  cv_I_pred <- (sd_new_subject / mu_new_subject) * 100
  mean_beta <- mean(fit$beta) 
  mean_df_I <- mean(fit$df_I)
  mean_df_A <- mean(fit$df_A)
  
  # Helper functions
  t_scale_to_sd <- function(scale, df) {
    scale * sqrt(df / (df - 2))
  }
  calculate_cv <- function(sd, mean) {
    (sd / mean) * 100
  }
  
  beta_pred <- fit$beta
  
  # These are population-average estimates
  cv_I <- calculate_cv(t_scale_to_sd(fit$sigma_I, mean_df_I), mean_beta)
  cv_A <- calculate_cv(t_scale_to_sd(fit$sigma_A, mean_df_A), mean_beta)
  cv_G <- calculate_cv(fit$sigma_G, mean_beta)
  
  # Calculate individual-specific CVs
  mu_i <- colMeans(fit$G) + mean_beta
  sds_i <- t_scale_to_sd(fit$sigma_i, mean_df_I)
  cv_pi <- t(t(sds_i) / mu_i) * 100
  
  return(list(
    beta_pred = beta_pred, 
    cv_I = cv_I, 
    cv_A = cv_A, 
    cv_G = cv_G,
    cv_I_pred = cv_I_pred,
    cv_pi = cv_pi
  ))
}

#' Process and Summarize Stan Model Output
#'
#' Takes the posterior draws from a Stan model fit and computes key biological
#' and analytical metrics, such as coefficients of variation (CVs) and credible
#' intervals, formatting them into user-friendly tables.
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
#' @param material A \code{character} string to label the biological material.
#' @param sex A \code{character} string to label the sex. Defaults to "Unknown".
#' @param group A \code{character} string used to label the group or cohort in
#'   the output tables. Defaults to "Y".
#' @param data An optional \code{data.frame} or \code{data.table} containing the
#'   original input data. If provided, its \code{SubjectID} column will be used to
#'   label the subjects in the subject-specific output table.
#'   Defaults to \code{NULL}.
#' @param digits An \code{integer} specifying the number of decimal places to use
#'   when formatting the credible intervals. Defaults to 1.
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
#' @section Calculation of CVs:
#' The method for calculating CVs depends on the \code{log_transformed} flag.
#' \itemize{
#'   \item \strong{If \code{log_transformed = TRUE}:} The function assumes that
#'     variance components (\code{sigma} parameters) are on the natural log scale.
#'     It uses the appropriate formulas to convert these into CVs on the original
#'     scale. For parameters with Student-t distributions (\code{sigma_I}, \code{sigma_A}),
#'     it accounts for the degrees of freedom (\code{df}).
#'   \item \strong{If \code{log_transformed = FALSE}:} The function assumes a model on the
#'     identity scale. CVs are calculated as \code{sd / mean * 100}. The standard
#'     deviation for t-distributed parameters is derived from the scale parameter
#'     (\code{sigma}) and degrees of freedom (\code{df}).
#' }
#'
#' @section Important Note on Model Validity:
#' The reliability of the metrics produced by this function depends entirely
#' on the quality and convergence of the Stan model fit. Before trusting the
#' output, it is essential to perform standard MCMC diagnostics (e.g., check for
#' divergent transitions, R-hat < 1.01, and sufficient \code{n_eff}).
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
process_stan_output <- function(fit, log_transformed, analyte = "X", material = "Unknown", sex = "Unknown",  group = "Y", data = NULL, digits = 1L) {
  
  # This central step keeps the complex math separate from table formatting.
  if (log_transformed) {
    cvs <- calculate_cvs_log_scale(fit)
  } else {
    cvs <- calculate_cvs_identity_scale(fit)
  }
  
  # --- Build Summary Table (Table 1) ---
  table1 <- data.table(
    "Measurand" = analyte,
    "Sex" = sex,          
    "Material" = material,
    "Group" = group,
    "No. of individuals" = ncol(fit$G_raw),
    "No. of results" = ncol(fit$I_raw),
    "Mean concentration (95% CrI)" = format_ci(c(mean(cvs$beta_pred), quantile(cvs$beta_pred, c(0.025, 0.975))), digits),
    "CVI (95% CrI) %" = format_ci(c(mean(cvs$cv_I), quantile(cvs$cv_I, c(0.025, 0.975))), digits),
    "dCV_P(i) 50% (20%-80%) %" = format_ci(quantile(cvs$cv_I_pred, c(0.50, 0.20, 0.80)), digits),
    "CVG (95% CrI) %" = format_ci(c(mean(cvs$cv_G), quantile(cvs$cv_G, c(0.025, 0.975))), digits),
    "CVA (95% CrI) %" = format_ci(c(mean(cvs$cv_A), quantile(cvs$cv_A, c(0.025, 0.975))), digits),
    "HBHR %" = format(sd(cvs$cv_I_pred) / mean(cvs$cv_I_pred) * 100, nsmall = 1, digits = 1)
  )
  
  # --- Build Subject-Specific Table (Table 2) ---
  individual_cvs_long <- data.table::melt(
    data.table(cvs$cv_pi),
    measure.vars = 1:ncol(cvs$cv_pi),
    variable.name = "SubjectIndex",
    value.name = "CV_pi"
  )
  # Convert V1, V2 -> 1, 2
  individual_cvs_long[, SubjectIndex := as.integer(gsub("V", "", SubjectIndex))]
  
  # Summarize by subject
  table2 <- individual_cvs_long[, .(
    `mean_CV_P(i)` = mean(CV_pi),
    `median_CV_P(i)` = median(CV_pi),
    `CV_P(i)_lwr` = quantile(CV_pi, 0.025),
    `CV_P(i)_upr` = quantile(CV_pi, 0.975)
  ), by = SubjectIndex]
  
  # Add population predictive distribution for plotting
  table2[, `:=`(
    `dCV_P(i)_median` = median(cvs$cv_I_pred),
    `dCV_P(i)_lwr` = quantile(cvs$cv_I_pred, 0.20),
    `dCV_P(i)_upr` = quantile(cvs$cv_I_pred, 0.80)
  )]
  
  # Add subject labels if provided
  if (!is.null(data)) {
    subject_ids <- unique(as.data.table(data)[, "SubjectID", with = FALSE])
    table2[, SubjectID := as.character(subject_ids[SubjectIndex, 1, with=FALSE][[1]])]
  } else {
    table2[, SubjectID := as.character(SubjectIndex)]
  }
  
  # Order by median CV and create an ordered factor for plotting
  data.table::setorder(table2, `median_CV_P(i)`)
  table2[, SubjectID := factor(SubjectID, levels = unique(SubjectID))]
  data.table::setcolorder(table2, c("SubjectIndex", "SubjectID"))
  
  
  # --- 4. Build Degrees of Freedom Table (Table 3) ---
  table3 <- data.table(
    "Measurand" = analyte,
    "Sex" = sex,          
    "Material" = material,
    "Group" = group,
    "No. of individuals" = ncol(fit$G_raw),
    "No. of results" = ncol(fit$I_raw),
    "Mean DFI (95% CrI)" = format_ci(c(mean(fit$df_I), quantile(fit$df_I, c(0.025, 0.975))), digits),
    "Mean DFA (95% CrI)" = format_ci(c(mean(fit$df_A), quantile(fit$df_A, c(0.025, 0.975))), digits)
  )
  
  # --- 5. Return Results ---
  return(list("Summary" = table1, "Subject_Specific" = table2, "Degrees_of_Freedom" = table3))
}

#' Plot Subject-Specific Intra-Individual CVs
#'
#' Creates a dot-and-whisker plot visualizing the 95% credible intervals for
#' the intra-individual CV (`CVI`) for each subject. It includes reference lines
#' for the population predictive distribution and allows for flexible grouping
#' by color and shape.
#'
#' @param processed_output A \code{list} object returned by
#'   \code{process_stan_output}. The function specifically uses the
#'   `Subject_Specific` data.table.
#' @param data An optional \code{data.frame} or \code{data.table} containing the
#'   original input data. This is **required** if you want to use the
#'   \code{color_by} or \code{shape_by} arguments, as it must contain the
#'   grouping columns.
#' @param color_by A \code{character} string specifying the column name in `data`
#'   to use for coloring the points (e.g., `"sex"`). Defaults to `NULL` (no color grouping).
#' @param shape_by A \code{character} string specifying the column name in `data`
#'   to use for shaping the points (e.g., `"material"`). Defaults to `NULL` (no shape grouping).
#' @param title A \code{character} string for the plot title. Defaults to a
#'   standard title.
#' @param subtitle A \code{character} string for the plot subtitle. Can be used
#'   to provide context like the analyte name.
#'
#' @details
#' The plot is designed for clear visual comparison of individual CVI estimates.
#' -   **Points and Error Bars:** The central point for each subject represents their
#'     median posterior CVI, and the horizontal bars represent the 95% credible interval.
#' -   **Green Reference Band:** The shaded green area represents the 60% prediction
#'     interval for a new subject's CVI (from the 20th to the 80th percentile).
#'     The solid green line is the median of this predictive distribution.
#' -   **Ordering:** Subjects are ordered from bottom to top by their median CVI.
#'
#' @return A \code{ggplot} object, which can be further customized.
#'
#' @export
plot_subject_specific_CVI <- function(processed_output, data = NULL,
                                      color_by = NULL, shape_by = NULL,
                                      title = "Subject-Specific CVs with 95% Credible Intervals",
                                      subtitle = NULL) {
  
  # --- 1. Input Validation and Data Preparation ---
  if (!"Subject_Specific" %in% names(processed_output)) {
    stop("Input 'processed_output' must be a list from the process_stan_output() function.")
  }
  
  plotting_data <- data.table::copy(processed_output$Subject_Specific)
  
  grouping_vars <- c(color_by, shape_by)
  if (length(grouping_vars) > 0) {
    if (is.null(data)) {
      stop("You must provide the 'data' argument to use 'color_by' or 'shape_by'.")
    }
    group_data <- as.data.table(data)
    missing_vars <- setdiff(grouping_vars, names(group_data))
    if (length(missing_vars) > 0) {
      stop(paste("Grouping variables not in 'data' frame:", paste(missing_vars, collapse = ", ")))
    }
    subject_metadata <- unique(group_data[, c("SubjectID", grouping_vars), with = FALSE])
    
    plotting_data[, SubjectID := as.character(SubjectID)]
    subject_metadata[, SubjectID := as.character(SubjectID)]
    
    plotting_data <- merge(plotting_data, subject_metadata, by = "SubjectID", all.x = TRUE)
  }
  
  # --- FIX: Re-sort the data by median CV AFTER the merge ---
  # This ensures the y-axis is ordered correctly from lowest to highest CV.
  data.table::setorder(plotting_data, `median_CV_P(i)`)
  
  # Create the ordered factor for the y-axis based on the new, correct sort order
  plotting_data[, SubjectID := factor(SubjectID, levels = unique(SubjectID))]
  
  # --- 2. Build the ggplot object ---
  p <- ggplot(data = plotting_data, aes(y = SubjectID))
  
  pop_summary <- unique(plotting_data[, .(`dCV_P(i)_lwr`, `dCV_P(i)_upr`, `dCV_P(i)_median`)])
  p <- p +
    geom_rect(
      aes(xmin = pop_summary$`dCV_P(i)_lwr`, xmax = pop_summary$`dCV_P(i)_upr`, ymin = -Inf, ymax = Inf),
      fill = "#28A745", alpha = 0.3, inherit.aes = FALSE
    ) +
    geom_vline(
      xintercept = c(pop_summary$`dCV_P(i)_lwr`,
                     pop_summary$`dCV_P(i)_median`,
                     pop_summary$`dCV_P(i)_upr`),
      color = "#1a702e", linewidth = 1.0
    )
  
  p <- p + geom_errorbarh(
    aes(xmin = `CV_P(i)_lwr`, xmax = `CV_P(i)_upr`),
    height = 0.5, linewidth = 0.75, color = "gray20"
  )
  
  # --- 3. Dynamically add point aesthetics ---
  point_aes <- aes(x = `median_CV_P(i)`)
  if (!is.null(color_by)) {
    point_aes$fill <- as.symbol(color_by)
  }
  if (!is.null(shape_by)) {
    point_aes$shape <- as.symbol(shape_by)
  }
  
  p <- p + geom_point(
    mapping = point_aes,
    size = 3,
    color = "gray20",
    shape = if (is.null(shape_by)) 21 else NA
  )
  
  # --- 4. Customize scales, labels, and theme ---
  if (!is.null(color_by)) {
    # Ensure the grouping variable is treated as a factor for discrete colors
    plotting_data[, (color_by) := as.factor(get(color_by))]
    p <- p + scale_fill_brewer(palette = "Set2", name = toTitleCase(color_by))
  } else {
    p <- p + scale_fill_manual(values = "orange", guide = "none")
  }
  
  if (!is.null(shape_by)) {
    plotting_data[, (shape_by) := as.factor(get(shape_by))]
    p <- p + scale_shape_manual(
      values = c(21, 22, 24, 25),
      name = toTitleCase(shape_by)
    )
  }
  
  p <- p +
    scale_x_continuous(
      name = expression(paste("Within-Individual CV (", CV[p(i)], ", %)")),
      labels = function(x) paste0(format(x, nsmall = 1, digits = 1), " %"),
      n.breaks = 8
    ) +
    ylab("Subject") +
    labs(title = title, subtitle = subtitle) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    )
  
  return(p)
}




