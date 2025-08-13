#include <Rcpp.h>
#include <cmath>      // For std::abs, sqrt, log, exp, std::isinf
#include <string>     // For std::string
#include <limits>     // For std::numeric_limits, needed for infinity()

// Helper function for lower-truncated normal distribution (simple rejection sampling)
// Truncates at zero
double rtnorm_lower(double mean, double sd, double lower_bound = 0.0) {
  if (sd < 0.0) {
    Rcpp::stop("SD must be non-negative in rtnorm_lower.");
  }
  // Check if sd is extremely small
  if (sd <= 1e-12) {
    return (mean >= lower_bound) ? mean : lower_bound;
  }
  
  double val;
  unsigned int tries = 0;
  const unsigned int max_tries = 10000; // Safety break
  
  do {
    Rcpp::checkUserInterrupt(); // Allow user to interrupt long C++ loops
    val = R::rnorm(mean, sd);
    if (++tries > max_tries) {
      Rcpp::warning("rtnorm_lower: Exceeded max_tries in rejection sampling. "
                      "Parameters (mean=%.2f, sd=%.2f, lower_bound=%.2f) "
                      "may lead to very low acceptance rates. "
                      "Returning boundary or last valid draw.", mean, sd, lower_bound);
      return (val < lower_bound) ? lower_bound : val; 
    }
  } while (val < lower_bound);
  return val;
}

//' Simulate Biological Variation Data with Extended Options
//'
//' @title Simulate Biological Variation Data From NTT Model
//' @name simulate_bv_data_ntt
//'
//' @param n \code{Integer}. The number of subjects.
//' @param S \code{Integer}. The number of samples taken from each subject. Defaults to 10.
//' @param R \code{Integer}. The number of replicates measured for each sample. Defaults to 2.
//' @param cvi \code{Double}. The coefficient of variation (in percent) for variability within each subject (default is 10).
//' @param cva \code{Double}. The coefficient of variation (in percent) for the variability of analytical error (default is 2).
//' @param cvg \code{Double}. The coefficient of variation (in percent) for variability between subjects (default is 50).
//' @param mu \code{Double}. The overall true mean for the population (default is 100). Must be non-negative.
//' @param hbhr \code{Double}. The Harris-Brown heterogeneity ratio (in percent) for the variability of within-subject standard deviations. Default is 0.
//' @param grand_mean \code{Boolean}. If \code{TRUE}, SD-components are based on `mu`. Otherwise, on subject/sample-specific means. Default is \code{TRUE}.
//' @param subject_dist_type \code{String}. Distribution for subject means: "normal" (results in folded normal via `std::abs()`), "truncnormal" (truncated at 0), or "lognormal". Default is "normal".
//' @param dfi \code{Double}. Degrees of freedom for Student's t-distribution of samples within subjects. Use `R_PosInf` (or `Inf`) for normal distribution. Default is `R_PosInf`. Must be > 0.
//' @param dfa \code{Double}. Degrees of freedom for Student's t-distribution of analytical replicates. Use `R_PosInf` (or `Inf`) for normal distribution. Default is `R_PosInf`. Must be > 0.
//'
//' @description
//' This function simulates biological variation data with options for different distributions
//' for subject means, within-subject variation, and analytical variation.
//' 
//' @details
//' \strong{Subject-Specific Means (`subject_means[i]`):}
//' Generated based on `mu`, `cvg`, and `subject_dist_type`:
//' \itemize{
//'   \item \code{"normal"}: `subject_means[i] = std::abs(R::rnorm(mu, sd_g))`, where `sd_g = (cvg/100) * mu`. This creates a folded normal distribution.
//'   \item \code{"truncnormal"}: `subject_means[i]` is drawn from `N(mu, sd_g)` truncated at 0. `sd_g = (cvg/100) * mu`. Values are inherently non-negative.
//'   \item \code{"lognormal"}: `subject_means[i]` is drawn from a log-normal distribution such that its arithmetic mean is `mu` and arithmetic CV is `cvg/100`. Values are inherently non-negative. `mu` must be > 0.
//' }
//' All `subject_means[i]` are thus non-negative.
//'
//' \strong{Samples within Subjects (`sample_val`):}
//' Drawn from a location-scale Student's t-distribution (or normal if `dfi = Inf`).
//' Location: `current_subject_mean`.
//' The "nominal SD" (`sdi_nominal`) is determined by `cvi` and potentially `hbhr`.
//'   - `mean_for_sdi_calc = grand_mean ? mu : current_subject_mean`
//'   - `mean_nominal_sdi = (cvi/100) * mean_for_sdi_calc`
//'   - `sd_for_nominal_sdi_draw = mean_nominal_sdi * (hbhr/100)`
//'   - `sdi_nominal = std::abs(R::rnorm(mean_nominal_sdi, sd_for_nominal_sdi_draw))`
//' Scale parameter for t-distribution (`actual_scale_s`):
//' \itemize{
//'   \item If `dfi == Inf` (normal): `actual_scale_s = sdi_nominal`. Draw from `N(current_subject_mean, actual_scale_s)`.
//'   \item If `dfi > 2`: `actual_scale_s = sdi_nominal / sqrt(dfi / (dfi - 2.0))`. Draw `current_subject_mean + actual_scale_s * R::rt(1, dfi)`. `sdi_nominal` is the target SD.
//'   \item If `0 < dfi <= 2`: `actual_scale_s = sdi_nominal`. Draw `current_subject_mean + actual_scale_s * R::rt(1, dfi)`. `sdi_nominal` is used as the scale parameter directly; the resulting distribution has infinite variance (or undefined for `dfi <= 1`).
//' }
//' Then, `sample_val = std::abs(raw_sample_val)`.
//'
//' \strong{Replicates within Samples (`replicate_val`):}
//' Drawn similarly using `cva`, `dfa`, and `sample_val` (or `mu` if `grand_mean`).
//' The "nominal SD" (`sva_nominal`) is `(cva/100) * (grand_mean ? mu : sample_val)`.
//' Scale parameter for t-distribution (`actual_scale_a`) is determined analogously to `actual_scale_s` using `sva_nominal` and `dfa`.
//' Then, `replicate_val = std::abs(raw_replicate_val)`.
//'
//' All final measured values (`y`) are non-negative due to `std::abs()`, resulting in folded-t or folded-normal distributions.
//'
//' @return A `DataFrame` with columns: `SubjectID`, `SampleID`, `ReplicateID`, `y`.
//'
//' @examples
//' \dontrun{
//' # Default (normal distributions for all components)
//' df_norm <- simulate_bv_data_extended(15, 10, 2, mu = 100)
//' summary(df_norm$y)
//'
//' # Student's t for samples (dfi=3) and replicates (dfa=5)
//' # R_PosInf can be passed as Inf from R
//' df_t <- simulate_bv_data_extended(15, 10, 2, mu = 100, dfi = 3, dfa = 5)
//' summary(df_t$y)
//'
//' # Log-normal subject means
//' df_lognorm_subj <- simulate_bv_data_extended(15, 10, 2, mu = 100,
//'                                            subject_dist_type = "lognormal")
//' summary(df_lognorm_subj$y)
//'
//' # Truncated normal subject means, t-dist for samples, grand_mean=FALSE
//' df_mix <- simulate_bv_data_extended(n = 50, S = 5, R = 2,
//'                                   cvi = 15, cva = 3, cvg = 60, mu = 50,
//'                                   hbhr = 20, grand_mean = FALSE,
//'                                   subject_dist_type = "truncnormal", dfi = 4.0, dfa = Inf)
//' summary(df_mix$y)
//' }
//'
// [[Rcpp::export]]
Rcpp::List simulate_bv_data_ntt(
   int n, int S = 10, int R = 2,
   double cvi = 10.0, double cva = 2.0, double cvg = 50.0,
   double mu = 100.0, double hbhr = 0.0, bool grand_mean = true,
   std::string subject_dist_type = "normal",
   double dfi = 99999999,
   double dfa = 99999999) {
 
 // --- Input Validation ---
 if (n <= 0 || S <= 1 || R <= 1) {
   Rcpp::stop("n must be positive, and S, and R must be at least 2L.");
 }
 if (mu < 0.0) { 
   Rcpp::stop("mu (overall true mean) must be non-negative.");
 }
 if (cvi < 0.0 || cva < 0.0 || cvg < 0.0 || hbhr < 0.0) {
   Rcpp::stop("CV and HBHR parameters must be non-negative.");
 }
 if (subject_dist_type != "normal" && subject_dist_type != "truncnormal" && subject_dist_type != "lognormal") {
   Rcpp::stop("subject_dist_type must be 'normal', 'truncnormal', or 'lognormal'.");
 }
 if (subject_dist_type == "lognormal" && mu <= 0.0) {
   Rcpp::stop("mu must be > 0 for lognormal distribution of subject means.");
 }
 if (dfi < 2.0 || dfa < 2.0) {
   Rcpp::stop("Degrees of freedom (dfi, dfa) must be at least 2.");
 }
 
 // --- Pre-calculations ---
 int total_size = n * S * R;
 Rcpp::IntegerVector SubjectID_vec(total_size);
 Rcpp::IntegerVector SampleID_vec(total_size);
 Rcpp::IntegerVector ReplicateID_vec(total_size);
 Rcpp::NumericVector y_vec(total_size);
 
 // --- CVs as decimal numbers ---
 double cvi_f = cvi / 100.0;
 double cva_f = cva / 100.0;
 double cvg_f = cvg / 100.0;
 double hbhr_f = hbhr / 100.0;
 
 // --- Generate subject means ---
 std::vector<double> subject_means(n);
 double sd_g = cvg_f * mu;
 
 if (subject_dist_type == "lognormal") {
   double var_arith_g = cvg_f * cvg_f;
   double log_sd_g = (cvg_f > 1e-9) ? std::sqrt(std::log(var_arith_g + 1.0)) : 0.0;
   double log_mean_g = std::log(mu) - 0.5 * log_sd_g * log_sd_g;
   for(int i = 0; i < n; i++) {
     subject_means[i] = std::exp(R::rnorm(log_mean_g, log_sd_g));
   }
 } 
 else if (subject_dist_type == "truncnormal") {
   for(int i = 0; i < n; i++) {
     subject_means[i] = rtnorm_lower(mu, sd_g, 0.0); // Ensures >= 0
   }
 } 
 else { // "normal" (implies folded normal via std::abs)
   for(int i = 0; i < n; i++) {
     subject_means[i] = std::abs(R::rnorm(mu, sd_g));
   }
 }
 
 // --- Generate samples and replicates ---
 int idx = 0;
 for(int i = 0; i < n; i++) { // Loop over subjects
   Rcpp::checkUserInterrupt();
   double current_subject_mean = subject_means[i];
   
   // Determine sdi_nominal (within-subject SD for samples)
   double sdi_mean_calc_base = grand_mean ? mu : current_subject_mean;
   double mean_nominal_sdi = cvi_f * sdi_mean_calc_base;
   double sd_for_sdi_draw = mean_nominal_sdi * hbhr_f;
   double sdi_nominal = std::abs(R::rnorm(mean_nominal_sdi, sd_for_sdi_draw));
   if (sdi_nominal < 1e-9 * mean_nominal_sdi && mean_nominal_sdi > 1e-9) { 
     sdi_nominal = 1e-9 * mean_nominal_sdi; 
   }
   if (sdi_nominal < 0) sdi_nominal = 0;
   
   // Determine actual scale for t-distribution for samples
   double actual_scale_s;
   if (dfi >= 99999999) { // Normal distribution
     actual_scale_s = sdi_nominal;
   } 
   else if (dfi > 2.0) {
     actual_scale_s = sdi_nominal / std::sqrt(dfi / (dfi - 2.0));
   } 
   else { // Will never happen
     actual_scale_s = sdi_nominal;
   }
   if (actual_scale_s < 0) actual_scale_s = 0;
   
   for(int j = 0; j < S; j++) { // Loop over samples
     double sample_val_raw;
     if (dfi >= 99999999) {
       sample_val_raw = R::rnorm(current_subject_mean, actual_scale_s);
     } 
     else {
       sample_val_raw = current_subject_mean + actual_scale_s * R::rt(dfi);
     }
     double sample_val = std::abs(sample_val_raw);
     
     // Determine sva_nominal (analytical SD for replicates)
     double sva_mean_calc_base = grand_mean ? mu : sample_val;
     double sva_nominal = cva_f * sva_mean_calc_base;
     if (sva_nominal < 0) sva_nominal = 0;
     
     // Determine actual scale for t-distribution for replicates
     double actual_scale_a;
     if (dfa >= 99999999) { // Normal distribution
       actual_scale_a = sva_nominal;
     } 
     else if (dfa > 2.0) { // t-dist with df > 2
       actual_scale_a = sva_nominal / std::sqrt(dfa / (dfa - 2.0));
     } 
     else { // Will never happen
       actual_scale_a = sva_nominal;
     }
     if (actual_scale_a < 0) actual_scale_a = 0;
     
     for(int k = 0; k < R; k++) { // Loop over replicates
       double replicate_val_raw;
       if (dfa >= 99999999) {
         replicate_val_raw = R::rnorm(sample_val, actual_scale_a);
       } 
       else {
         replicate_val_raw = sample_val + actual_scale_a * R::rt(dfa);
       }
       double replicate_val = std::abs(replicate_val_raw);
       
       SubjectID_vec[idx] = i + 1;
       SampleID_vec[idx] = j + 1;
       ReplicateID_vec[idx] = k + 1;
       y_vec[idx] = replicate_val;
       idx++;
     }
   }
 }
 
 return Rcpp::List::create(
   Rcpp::Named("SubjectID") = SubjectID_vec,
   Rcpp::Named("SampleID") = SampleID_vec,
   Rcpp::Named("ReplicateID") = ReplicateID_vec,
   Rcpp::Named("y") = y_vec);
}