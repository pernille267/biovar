#include <Rcpp.h>
using namespace Rcpp;


//' @title Generate Stan Data List for Priors
//' 
//' @name process_stan_data_priors
//'
//' @param beta A \code{double} representing the expected value of the main parameter.
//' @param cvi A \code{double} for the expected Coefficient of Variation for
//'   individual-level variance (I), provided as a percentage (e.g., 10 for 10%).
//'   Must be strictly positive. Defaults to 10.
//' @param cva A \code{double} for the expected Coefficient of Variation for
//'   assay-level variance (A), provided as a percentage. Must be strictly
//'   positive. Defaults to 3.
//' @param cvg A \code{double} for the expected Coefficient of Variation for
//'   group-level variance (G), provided as a percentage. Must be strictly
//'   positive. Defaults to 20.
//' @param dfi A \code{double} for the expected degrees of freedom for the
//'   individual-level t-distribution. Must be >= 2.1. Defaults to 9999.
//' @param dfa A \code{double} for the expected degrees of freedom for the
//'   assay-level t-distribution. Must be >= 2.1. Defaults to 9999.
//' @param strength A \code{NumericVector} of length 6 containing non-negative
//'   multipliers that control the "tightness" (standard deviation) of the
//'   hyperpriors. The elements correspond to the priors for:
//'   (1) \code{beta}, 
//'   (2) \code{sigma_i} (mean), 
//'   (3) \code{sigma_A},
//'   (4) \code{sigma_G},
//'   (5) \code{df_I} and
//'   (6) \code{df_A}.
//'   Defaults to \code{c(1, 1, 1, 1, 1, 1)}.
//' @param log_transformed A \code{bool} indicating whether the underlying model
//'   assumes the data is log-transformed. If \code{TRUE}, priors for standard
//'   deviations (\code{sigma_I}, \code{sigma_A}, \code{sigma_G}) are calculated on the log scale.
//'   Defaults to \code{FALSE}.
//' 
//' @description
//' This function processes user-friendly prior specifications and converts them
//' into a list of hyperparameters (means and standard deviations) suitable for
//' a Stan model.
//' 
//' @details
//' This function processes user-friendly prior specifications and converts them
//' into a list of hyperparameters (means and standard deviations) suitable for
//' a Stan model. It takes an expected value (\code{beta}), coefficients of variation
//' (CVs), and degrees of freedom, and calculates the corresponding parameters
//' for the prior distributions.
//' 
//' @return A \code{list} containing the calculated hyperparameters for the Stan model's priors:
//' \describe{
//'   \item{\code{prior_beta_mean}}{Mean for the prior on \code{beta}.}
//'   \item{\code{prior_beta_sd}}{Standard deviation for the prior on \code{beta}.}
//'   \item{\code{prior_sigma_i_mean_mean}}{Mean for the hyperprior on the mean of \code{sigma_i}.}
//'   \item{\code{prior_sigma_i_mean_sd}}{SD for the hyperprior on the mean of \code{sigma_i}.}
//'   \item{\code{prior_sigma_i_sd_mean}}{Mean for the hyperprior on the SD of \code{sigma_i}.}
//'   \item{\code{prior_sigma_i_sd_sd}}{SD for the hyperprior on the SD of \code{sigma_i}.}
//'   \item{\code{prior_sigma_A_mean}}{Mean for the prior on \code{sigma_A}.}
//'   \item{\code{prior_sigma_A_sd}}{Standard deviation for the prior on \code{sigma_A}.}
//'   \item{\code{prior_sigma_G_mean}}{Mean for the prior on \code{sigma_G}.}
//'   \item{\code{prior_sigma_G_sd}}{Standard deviation for the prior on \code{sigma_G}.}
//'   \item{\code{prior_df_I_mean}}{Mean for the prior on \code{df_I}.}
//'   \item{\code{prior_df_I_sd}}{Standard deviation for the prior on \code{df_I}.}
//'   \item{\code{prior_df_A_mean}}{Mean for the prior on \code{df_A}.}
//'   \item{\code{prior_df_A_sd}}{Standard deviation for the prior on \code{df_A}.}
//' }
//'
//' @examples
//' # Generate priors with default settings for a beta of 100
//' process_stan_data_priors(beta = 100)
//'
//' # Generate priors for a log-transformed model with tighter CVs
//' process_stan_data_priors(beta = 7, cvi = 5, cva = 2, log_transformed = TRUE)
//'
//' # Increase the uncertainty (standard deviation) of the beta prior
//' process_stan_data_priors(beta = 100, strength = c(2.0, 1, 1, 1, 1, 1))
//'

// [[Rcpp::export]]
List process_stan_data_priors(double beta, double cvi = 10, double cva = 3, double cvg = 20, double dfi = 9999, double dfa = 9999, NumericVector strength = NumericVector::create(1, 1, 1, 1, 1, 1), bool log_transformed = false) {
  int strength_length = strength.size();
  
  // Validity checks
  if (strength_length != 6) {
    Rcpp::stop("strength must be of length 6.");
  }
  if (cvi <= 0) {
    Rcpp::stop("cvi (hyper for E[CV_I]) must be strictly positive");
  }
  if (cva <= 0) {
    Rcpp::stop("cva (hyper for E[CV_A]) must be strictly positive");
  }
  if (cvg <= 0) {
    Rcpp::stop("cvg (hyper for E[CV_G]) must be strictly positive");
  }
  if (dfi <= 2.1) {
    Rcpp::stop("dfi (hyper for E[df_I]) must be larger than or equal to 2.1");
  }
  if (dfa <= 2.1) {
    Rcpp::stop("dfa (hyper for E[df_A]) must be larger than or equal to 2.1");
  }
  
  for (int i = 1; i < strength_length; ++i) {
    if (strength[i] < 0) {
      Rcpp::stop("Each element of strength must be non-negative");
    }
  }
  
  // CV as decimal numbers
  double CV_I_prior = cvi / 100.0;
  double CV_A_prior = cva / 100.0;
  double CV_G_prior = cvg / 100.0;
  
  double beta_prior = 0;
  double sigma_I_prior = 0;
  double sigma_A_prior = 0;
  double sigma_G_prior = 0;
  
  if (log_transformed) {
    beta_prior = std::log(beta);
    sigma_I_prior = sqrt(std::log(pow(CV_I_prior, 2.0) + 1));
    sigma_A_prior = sqrt(std::log(pow(CV_A_prior, 2.0) + 1));
    sigma_G_prior = sqrt(std::log(pow(CV_G_prior, 2.0) + 1));
  }
  else {
    beta_prior = beta;
    sigma_I_prior = CV_I_prior * beta_prior;
    sigma_A_prior = CV_A_prior * beta_prior;
    sigma_G_prior = CV_G_prior * beta_prior;
  }
  
  if (beta_prior < 0) {
    strength[0] = strength[0] * (-1.0);
  }
  
  return List::create(
    Named("prior_beta_mean") = beta_prior,
    Named("prior_beta_sd") = strength[0] * beta_prior,
    Named("prior_sigma_i_mean_mean") = sigma_I_prior,
    Named("prior_sigma_i_mean_sd") = strength[1] * sigma_I_prior,
    Named("prior_sigma_i_sd_mean") = 0.5 * sigma_I_prior,
    Named("prior_sigma_i_sd_sd") = 2 * sigma_I_prior,
    Named("prior_sigma_A_mean") = sigma_A_prior,
    Named("prior_sigma_A_sd") = strength[2] * sigma_A_prior,
    Named("prior_sigma_G_mean") = sigma_G_prior,
    Named("prior_sigma_G_sd") = strength[3] * sigma_G_prior,
    Named("prior_df_I_mean") = dfi,
    Named("prior_df_I_sd") = strength[4] * dfi,
    Named("prior_df_A_mean") = dfa,
    Named("prior_df_A_sd") = strength[5] * dfa
  );
  
}


