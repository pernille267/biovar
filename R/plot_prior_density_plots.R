#' Draw Prior Density Plots for Bayesian Models
#'
#' @param beta A \code{double}. The expected value of the fixed effect
#'             parameter, \eqn{\beta}.
#' @param cvi A \code{double}. The expected value of the mean Coefficient
#'            of Variation for the biological variation within subjects
#'            (\eqn{\mathrm{E}[\mathrm{E}[\mathrm{CV}_i]]}), provided as a
#'            percentage (e.g., 10 for 10%). Must be strictly positive.
#'            Defaults to 10.
#' @param cva A \code{double}. The expected Coefficient of Variation for
#'            IVD-MD measurement error
#'            (\eqn{\mathrm{E}[\mathrm{CV}_{\mathrm{A}}]}), provided as a
#'            percentage. Must be strictly positive. Defaults to 3.
#' @param cvg A \code{double}. The expected Coefficient of Variation for
#'            for the biological variation between subjects
#'            (\eqn{\mathrm{E}[\mathrm{CV}_{\mathrm{A}}]}), provided as a
#'            percentage. Must be strictly positive. Defaults to 20.
#' @param dfi A \code{double}. The expected degrees of freedom
#'            \eqn{\mathrm{E}[\mathrm{df}_{\mathrm{I}}]} for the within
#'            subjects effect
#'            (\eqn{\mathrm{I}_{is} \sim \mathrm{lst}(\mathrm{df}_{\mathrm{I}},
#'            0, \sigma_i^2)})
#'            of the assumed model. Must be \eqn{\geq 2}. Defaults to 20.
#' @param dfa A \code{double}. The expected degrees of freedom
#'            \eqn{\mathrm{E}[\mathrm{df}_{\mathrm{A}}]} for the IVD-MD
#'            measurement error effect
#'            (\eqn{\mathrm{A}_{isr} \sim \mathrm{lst}(\mathrm{df}_{\mathrm{A}},
#'            0, \sigma_A^2)})
#'            of the assumed model. Must be \eqn{\geq 2}. Defaults to 20.
#' @param strength A \code{numeric} vector of length \code{6}. These are
#'                 non-negative multipliers that control how informative the
#'                 priors of \code{beta}, \code{cvi}, \code{cva}, \code{cvg},
#'                 \code{dfi}, and \code{dfa} are required to be. For less
#'                 informative priors, choose larger values for the 'strength'
#'                 vector, and vice versa. The first element of strength
#'                 applies to \code{beta}, the second to \code{cvi} and so on
#'                 adhering to the order of parameters given above. Defaults to
#'                 unit strength, \code{rep(1, 6)}.
#' @param log_transformed A non-missing \code{logical} value. If \code{TRUE},
#'                        the multiplicative model,
#'                        \eqn{y_{isr} = \beta\cdot G_i \cdot I_{is} \cdot {A_{isr}}}
#'                        is assumed. In practice, this means that we model
#'                        \eqn{\log(y_{isr})}. The priors for \eqn{\beta} and
#'                        the variance components (\eqn{\sigma_i^2},
#'                        \eqn{\sigma_A}, and \eqn{\sigma_G}) are calculated on
#'                        the log scale. If \code{FALSE}, the standard additive
#'                        model is assumed. Defaults to \code{FALSE}.
#' @param model A \code{character} string. The name of the Bayesian model used.
#'              Currently, there are only two options, \code{NTT} and
#'              \code{NTTDFGAM}. The latter uses the gamma distribution a
#'              prior for \code{dfi} and \code{dfa}.
#'              
#' @description
#' Constructs and prints prior density plots for the parameters of the Bayesian
#' model.
#' 
#' @details
#' This function visualizes the prior distributions specified by the user's
#' input. It operates in several steps:
#' \enumerate{
#'    \item It calls \code{process_stan_data_priors()} to translate the
#'    user-friendly inputs (e.g., \code{cvi = 10} for 10\%) into
#'    hyperparameters (e.g., mean and standard deviation) of the 'effective'
#'    prior distributions (e.g., Truncated Normal, Gamma) that are actually
#'    used in the Stan model.
#'    \item It draws a large number of samples (\code{N = 1e5}) from each of
#'    these effective prior distributions.
#'    \item For hierarchical parameters like \code{cvi}, it samples from the
#'    full hierarchical structure 
#'    (e.g., \code{cvi} ~ TN(\code{cvi_mean}, \code{cvi_sd}), where
#'    \code{cvi_mean} and \code{cvi_sd} are themselves drawn from priors).
#'    \item It transforms these samples from the 'effective' scale back to the
#'    'user' scale (e.g., Coeffient of Variation in \%)
#'    \item It generates a \code{ggplot2} density plot for each parameter.
#'    \item It combines all 9 plots (for \code{beta}, \code{dfi}, \code{dfa},
#'    \code{cvi}, \code{cva}, \code{cvg}, \code{cvi_mean}, \code{cvi_sd}, and 
#'    \code{hbhr}) into a single plot using \code{cowplot::plot_grid}
#' }
#'
#' @returns
#' A \code{cowplot} object (which is also a \code{ggplot} object) containing a
#' 3 X 3 grid of all 9 prior density plots. When printed, this object will
#' display the combined plot.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # NOTE: These examples require the 'process_stan_data_priors' function
#' # to be exported and available in your package's namespace.
#'
#' # Example 1: Basic plot with default settings, specifying only beta
#' # Assumes an additive model (log_transformed = FALSE)
#' plot_priors_default <- plot_prior_density_plots(beta = 100)
#' print(plot_priors_default)
#'
#'
#' # Example 2: Log-transformed model with custom CVs and degrees of freedom
#' # Uses the "NTTDFGAM" model for priors on df
#' plot_priors_log_custom <- plot_prior_density_plots(
#'   beta = 6.0,            # E.g., log(400)
#'   cvi = 15,
#'   cva = 5,
#'   cvg = 25,
#'   dfi = 10,
#'   dfa = 10,
#'   log_transformed = TRUE,
#'   model = "NTTDFGAM"
#' )
#' print(plot_priors_log_custom)
#'
#'
#' # Example 3: Additive model with very uninformative (high strength) priors
#' plot_priors_uninformative <- plot_prior_density_plots(
#'   beta = 100,
#'   strength = rep(10, 6)
#' )
#' print(plot_priors_uninformative)
#' }
plot_prior_density_plots <- function(beta, cvi = 10, cva = 3, cvg = 20, dfi = 20, dfa = 20, strength = c(1, 1, 1, 1, 1, 1), log_transformed = FALSE, model = "NTT") {
  
  # --- Constants --------------------------------------------------------------
  N <- 1e5
  hbhr_mean <- 1/2
  hbhr_sd <- 2/3
  # ... Maybe we will add more constants (if necessary)
  
  # --- Helper functions -------------------------------------------------------
  
  # Draw from effective priors (those priors used in the stan scripts)
  sample_from_effective_prior <- function(N = 1e5, hyperparameters, prior = "TN", additional_prior_parameters = c(0, 0, Inf)) {
    
    # Track the samples
    samples_from_effective_prior <- NULL
    
    # Validation check of `hyperparameters` (for debugging mainly)
    # Check if hyperparameters is a list. If so, convert it to a vector
    if (is.list(hyperparameters)) {
      if (all(c("mean", "sd") %in% names(hyperparameters))) {
        hyperparameters <- c(hyperparameters$mean, hyperparameters$sd)  
      }
      else {
        if (length(hyperparameters) >= 2) {
          if (is.numeric(hyperparameters[[1]]) & is.numeric(hyperparameters[[2]])) {
            hyperparameters <- c(hyperparameters[[1]], hyperparameters[[2]])  
          }
          else {
            stop(
              "hyperparameters is a list with more than two elements.",
              "but the at least one of the first two elements are not numeric.",
              "This should not happen at all. This is your fault, Pernille."
            )
          }
        }
      }
      
    }
    # Throw an error if `hyperparameters` is not numeric
    if (!is.numeric(hyperparameters)) {
      stop(
        "hyperparameters must be a numeric vector of length two.",
        "This is your fault, Pernille. Not the user. Fix it."
      )
    }
    # Throw an error if `hyperparameters` is not a vector of length two
    if (length(hyperparameters) != 2) {
      stop(
        "hyperparameters must be of length two.",
        "This is your fault, Pernille. Not the user. Fix it."
      )
    }
    # Throw an error if the second element of hyperparameters is negative
    if (hyperparameters[2] < 0) {
      stop(
        "The standard deviation part of hyperparameters is negative.",
        "This is your fault Pernille. Fix it."
      )
    }
    
    # prior = TN: Sample from a truncated normal distribution
    # additional_prior_parameters: 
    # requirement: Length 3L
    # additional_prior_parameters[1]: constant shift to all generated samples
    # additional_prior_parameters[2]: Lower boundary for the truncated normal distribution (can be -Inf)
    # additional_prior_parameters[3]: Upper boundary for the truncated normal distribution (can be Inf)
    if (prior == "TN" | prior == "tn") {
      samples_from_effective_prior <- truncnorm::rtruncnorm(
        n = N,
        a = additional_prior_parameters[2],
        b = additional_prior_parameters[3],
        mean = hyperparameters[1],
        sd = hyperparameters[2]
      )
      samples_from_effective_prior <- samples_from_effective_prior + additional_prior_parameters[1]
    }
    # prior = N: Sample from a normal distribution
    # additional_prior_parameters: 
    # requirement: Length 1L
    # additional_prior_parameters[1]: constant shift to all generated samples
    else if (prior == "N" | prior == "n" | prior == "norm") {
      samples_from_effective_prior <- rnorm(
        n = N,
        mean = hyperparameters[1],
        sd = hyperparameters[2]
      )
      samples_from_effective_prior <- samples_from_effective_prior + additional_prior_parameters[1]
    }
    # prior = G: Sample from a gamma distribution
    # additional_prior_parameters: 
    # requirement: Length 1L
    # additional_prior_parameters[1]: constant shift to all generated samples
    else if (prior == "G" | prior == "g" | prior == "gamma") {
      g_alpha <- hyperparameters[1]^2 / hyperparameters[2]^2
      g_lambda <- hyperparameters[1] / hyperparameters[2]^2
      samples_from_effective_prior <- rgamma(
        n = N,
        shape = g_alpha,
        rate = g_lambda
      )
      samples_from_effective_prior <- samples_from_effective_prior + additional_prior_parameters[1]
    }
    
    if (is.null(samples_from_effective_prior)) {
      stop(
        "No samples were performed. Likely due to problems with the prior",
        "specification."
      )
    }
    
    return(samples_from_effective_prior)
    
    
    
  }
  
  # Transform to user-friendly priors (those priors that users work with)
  transform_effective_prior_samples_to_user_prior_samples <- function(samples, parameter = "beta", log_transformed = FALSE, additional_parameters = NULL) {
    # If parameter = beta, we can assume that the effective samples are drawn
    # from a normal distribution.
    if (parameter == "beta") {
      # If log(beta) ~ N(...,...)
      # We take the exponential: exp(log(beta)) = beta ~ lognormal(..., ...)
      if (log_transformed) {
        user_prior_samples <- exp(samples)
      }
      else {
        user_prior_samples <- samples
      }
    }
    # If parameter is `cvi_mean`, `cvi_sd`, `cvi`, `cva` or `cvg`
    # we can assume that the effective samples are drawn from a truncated
    # normal distribution.
    else if (any(parameter == c("cvi_mean", "cvi_sd", "cvi", "cva", "cvg"))) {
      # Extract information from additional_parameters
      # additional_parameters[1]: E[df_I] or E[df_A] (only if parameter == "cva"). Set to large number if parameter == "cvg".
      # additional_parameters[2]: E[beta] (mandatory, but only relevant if log_transformed == FALSE)
      lst_scale <- additional_parameters[1] / (additional_parameters[1] - 2)
      beta_scale <- additional_parameters[2]
      # If on the log-scale:
      # LSD = E[LSD_i], SD[LSD_i], LSD_i, LSD_A, LSDG
      # LSD ~ TN(mu_log, sigma_log, 0, Inf)
      # We transform from LSD (samples) to
      # CV = E[CV_i], SD[CV_i], CV_i, CV_A, CVG 
      # Using: CV = sqrt(exp(LSD^2 * E[df_*] / (E[df_*] - 2)) - 1) * 100
      if (log_transformed) {
        user_prior_samples <- sqrt(exp(samples^2 * lst_scale) - 1) * 100
      }
      # If on the identity-scale:
      # SD = E[LSD_i], SD[LSD_i], LSD_i, LSD_A, LSDG
      # SD ~ TN(mu, sigma, 0, Inf)
      # We transform from SD (samples) to
      # CV = E[CV_i], SD[CV_i], CV_i, CV_A, CVG 
      # Using: CV = SD * sqrt(E[df_*] / (E[df_*] - 2)) / E[beta] * 100
      else {
        user_prior_samples <- samples * sqrt(lst_scale) / beta_scale * 100
      }
    }
    else {
      user_prior_samples <- samples
    }
    
    return(user_prior_samples)
    
  }
  
  # Draw a density plot for a particular parameter
  draw_prior_density_plot <- function(user_samples, kernel_multiplier = 2, density_limits = c(0, Inf), vline_x = 0, soft_trim = 0.01, hard_trim = 0.0025, ticks = 8, xname = "x", attach_to_x_breaks_labels = "") {
    
    # --- Inferred plotting parameters -----------------------------------------
    
    # ------ Limits for x-values ---
    
    # --------- Soft limits ---
    soft_xlims <- quantile(user_samples, probs = c(soft_trim, 1 - soft_trim), names = FALSE)
    if (sum(user_samples < 2) == 0) { # Indicates that samples are for dfi or dfa (hopefully)
      soft_xlims[1] <- 2  
    }
    else if (sum(user_samples < 0) == 0) { # Indicates that samples are for variability components (hopefully)
      soft_xlims[1] <- 0
    }
    # --------- Calculation limits ---
    calc_xlims <- c(-NA, NA)
    if (soft_xlims[1] == 2) {
      calc_xlims[1] <- 2
    }
    else if (soft_xlims[1] == 0) {
      calc_xlims[1] <- 0
    }
    
    # --------- Hard limits ---
    hard_xlims <- c(-Inf, Inf)
    if (soft_xlims[1] < 0) {
      hard_xlims <- quantile(user_samples, probs = c(hard_trim, 1 - hard_trim), names = FALSE)
    }
    else {
      hard_xlims <- c(
        soft_xlims[1],
        quantile(user_samples, probs = 1 - hard_trim * 2, names = FALSE)
      )
    }
    
    # ------ Breaks on the x-axis ---
    x_breaks <- seq(
      from = soft_xlims[1],
      to = soft_xlims[2],
      length.out = ticks
    )
    
    # --------- Labels for the breaks on the x-axis ---
    x_breaks_labels <- x_breaks
    if (diff(range(x_breaks)) < ticks / 10) {
      x_breaks_labels <- format(
        x = x_breaks,
        nsmall = 2,
        digits = 2
      )
    }
    else if (diff(range(x_breaks)) < ticks) {
      x_breaks_labels <- format(
        x = x_breaks,
        nsmall = 1,
        digits = 1
      )
    }
    else {
      x_breaks_labels <- as.character(
        x = round(x_breaks)
      )
    }
    x_breaks_labels <- paste0(x_breaks_labels, attach_to_x_breaks_labels)
    
    # Construct density plot
    prior_density_plot <- ggplot() +
      geom_density(
        mapping = aes(
          x = user_samples
        ),
        adjust = kernel_multiplier,
        bounds = density_limits,
        color = "#605CA8",
        fill = "#28A745",
        alpha = 0.7
      ) +
      geom_vline(
        xintercept = vline_x
      ) +
      scale_x_continuous(
        name = xname,
        breaks = x_breaks,
        labels = x_breaks_labels,
        limits = calc_xlims,
        expand = expansion(mult = c(0, 0))
      ) +
      scale_y_continuous(
        name = "Density",
        breaks = NULL,
        expand = expansion(mult = c(0, 0.02))
      ) +
      coord_cartesian(
        xlim = hard_xlims
      ) +
      theme_classic()
  }
  
  # Calculate means and standard deviations for the seven prior distributions
  effective_prior_hyperparameters <- process_stan_data_priors(
    beta = beta,
    cvi = cvi * (dfi - 2) / dfi,
    cva = cva * (dfa - 2) / dfa,
    cvg = cvg,
    dfi = dfi,
    dfa = dfa,
    strength = strength,
    log_transformed = log_transformed
  )
  
  # Make list of hyperparameters
  hyperparameters_list <- list(
    "beta" = beta,
    "cvi_mean" = cvi,
    "cvi_sd" = hbhr_mean * cvi,
    "hbhr" = hbhr_mean * 100,
    "cvi" = cvi,
    "cva" = cva,
    "cvg" = cvg,
    "dfi" = dfi,
    "dfa" = dfa
  )
  
  # Make list of effective hyperparameters
  effective_hyperparameters_list <- list(
    "beta" = c(effective_prior_hyperparameters[[1]], effective_prior_hyperparameters[[2]]),
    "cvi_mean" = c(effective_prior_hyperparameters[[3]], effective_prior_hyperparameters[[4]]),
    "cvi_sd" = c(effective_prior_hyperparameters[[5]], effective_prior_hyperparameters[[6]]),
    "hbhr" = c(hbhr_mean, hbhr_sd) * 100,
    "cvi" = NA,
    "cva" = c(effective_prior_hyperparameters[[7]], effective_prior_hyperparameters[[8]]),
    "cvg" = c(effective_prior_hyperparameters[[9]], effective_prior_hyperparameters[[10]]),
    "dfi" = c(effective_prior_hyperparameters[[11]], effective_prior_hyperparameters[[12]]),
    "dfa" = c(effective_prior_hyperparameters[[13]], effective_prior_hyperparameters[[14]])
  )
  
  # Make list of prior distributions
  prior_distribution_list <- list(
    "beta" = "N",
    "cvi_mean" = "TN",
    "cvi_sd" = "TN",
    "hbhr" = "TN",
    "cvi" = "TN",
    "cva" = "TN",
    "cvg" = "TN",
    "dfi" = ifelse(model == "NTTDFGAM", "G", "TN"),
    "dfa" = ifelse(model == "NTTDFGAM", "G", "TN")
  )
  
  # Make list of additional prior distribution parameters
  additional_prior_parameters_list <- list(
    "beta" = c(0, -Inf, Inf),
    "cvi_mean" = c(0, 0, Inf),
    "cvi_sd" = c(0, 0, Inf),
    "hbhr" = c(0, 0, Inf),
    "cvi" = c(0, 0, Inf),
    "cva" = c(0, 0, Inf),
    "cvg" = c(0, 0, Inf),
    "dfi" = c(2, 0, Inf),
    "dfa" = c(2, 0, Inf)
  )
  
  # Make list of additional parameters
  additional_parameters_list <- list(
    "beta" = 0,
    "cvi_mean" = c(dfi, beta),
    "cvi_sd" = c(dfi, beta),
    "hbhr" = 0,
    "cvi" = c(dfi, beta),
    "cva" = c(dfa, beta),
    "cvg" = c(1e6, beta),
    "dfi" = 0,
    "dfa" = 0
  )
  
  # Sample from all user priors
  all_user_priors_samples <- sapply(
    X = 1:9,
    FUN = function(par_id) {
      if (any(is.na(effective_hyperparameters_list[[par_id]]))) {
        par_effective_prior_samples <- truncnorm::rtruncnorm(
          n = N,
          a = 0,
          b = Inf,
          mean = sample_from_effective_prior(
            N = N,
            hyperparameters = effective_hyperparameters_list[[2]],
            prior = "TN",
            additional_prior_parameters = additional_prior_parameters_list[[2]]
          ),
          sd = sample_from_effective_prior(
            N = N,
            hyperparameters = effective_hyperparameters_list[[3]],
            prior = "TN",
            additional_prior_parameters = additional_prior_parameters_list[[3]]
          )
        )
        par_user_prior_samples <- transform_effective_prior_samples_to_user_prior_samples(
          samples = par_effective_prior_samples,
          parameter = names(effective_hyperparameters_list)[par_id],
          log_transformed = log_transformed,
          additional_parameters = additional_parameters_list[[par_id]]
        )
      }
      else {
        par_user_prior_samples <- transform_effective_prior_samples_to_user_prior_samples(
          samples = sample_from_effective_prior(
            N = N,
            hyperparameters = effective_hyperparameters_list[[par_id]],
            prior = prior_distribution_list[[par_id]],
            additional_prior_parameters = additional_prior_parameters_list[[par_id]]
          ),
          parameter = names(effective_hyperparameters_list)[par_id],
          log_transformed = log_transformed,
          additional_parameters = additional_parameters_list[[par_id]]
        )
      }
      return(par_user_prior_samples)
    },
    simplify = FALSE
  )
  
  names(all_user_priors_samples) <- names(effective_hyperparameters_list)
  
  # Create one seperate density plot for each prior
  
  # Make list of density limits for each parameter
  density_limits_list <- list(
    "beta" = c(-Inf, Inf),
    "cvi_mean" = c(0, Inf),
    "cvi_sd" = c(0, Inf),
    "hbhr" = c(0, Inf),
    "cvi" = c(0, Inf),
    "cva" = c(0, Inf),
    "cvg" = c(0, Inf),
    "dfi" = c(2, Inf),
    "dfa" = c(2, Inf)
  )
  
  # Make list of vline x-coordinates
  vline_x_list <- sapply(
    X = 1:9,
    FUN = function(par_id) {
      effective_mean <- mean(all_user_priors_samples[[par_id]])
      user_mean <- hyperparameters_list[[par_id]]
      lower_bound <- c(0, 0, 0, 0, 0, 0, 0, 2, 2)[par_id]
      return(c(lower_bound, user_mean, effective_mean))
    },
    simplify = FALSE
  )
  
  # Make list of xnames
  xname_list <- list(
    "beta" = expression(beta),
    "cvi_mean" = expression(E(CV[i])),
    "cvi_sd" = expression(SD(CV[i])),
    "hbhr" = "HBHR",
    "cvi" = expression(CV[i]),
    "cva" = expression(CV[A]),
    "cvg" = expression(CV[G]),
    "dfi" = expression(df[I]),
    "dfa" = expression(df[A])
  )
  
  # Make list of break labels appendices (e.g., `%`)
  attach_to_x_breaks_labels_list <- list(
    "beta" = "",
    "cvi_mean" = " %",
    "cvi_sd" = " %",
    "hbhr" = " %",
    "cvi" = " %",
    "cva" = " %",
    "cvg" = " %",
    "dfi" = "",
    "dfa" = ""
  )
  
  # Construct the plots for each parameter
  prior_density_plots <- sapply(
    X = 1:9,
    FUN = function(par_id) {
      draw_prior_density_plot(
        user_samples = all_user_priors_samples[[par_id]],
        density_limits = density_limits_list[[par_id]],
        vline_x = vline_x_list[[par_id]],
        xname = xname_list[[par_id]],
        attach_to_x_breaks_labels = attach_to_x_breaks_labels_list[[par_id]] 
      )
    },
    simplify = FALSE
  )
  
  # Name the plots list for easier access
  names(prior_density_plots) <- names(effective_hyperparameters_list)
  
  # Gather all plots into one object using cowplot::plot_grid
  gathered_prior_plots <- cowplot::plot_grid(
    plotlist = list(
      prior_density_plots$beta,
      prior_density_plots$dfi,
      prior_density_plots$dfa,
      prior_density_plots$cvi,
      prior_density_plots$cva,
      prior_density_plots$cvg,
      prior_density_plots$cvi_mean,
      prior_density_plots$cvi_sd,
      prior_density_plots$hbhr
    )
  )
  
  gathered_prior_plots
  
}