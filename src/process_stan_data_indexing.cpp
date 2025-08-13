#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>
#include <numeric> // For std::accumulate
#include <algorithm> // For std::max

using namespace Rcpp;

//' Process Stan Data Indexing
//'
//' This function processes subject and sample data to create the necessary indexing
//' structures for hierarchical Bayesian models in Stan. It converts character-based
//' subject IDs to numeric indices and creates comprehensive mapping structures for
//' efficient Stan model execution.
//'
//' @param SubjectID_orig_R A \code{character} vector containing subject identifiers.
//'   Each element represents the subject ID for the corresponding observation.
//' @param SampleID_orig_R An \code{integer} vector containing sample identifiers within
//'   each subject. Must be positive \code{integer}s (1-based indexing).
//' @param y_R A numeric vector containing the response variable observations.
//'   This is passed through unchanged but must match the length of other inputs.
//'
//' @return A named list containing the following components for Stan model use:
//' \describe{
//'   \item{N_obs}{\code{integer}: Total number of observations}
//'   \item{y}{Numeric vector: Response variable (passed through unchanged)}
//'   \item{N_subj}{\code{integer}: Total number of unique subjects}
//'   \item{subj_idx}{\code{integer} vector: Subject index for each observation (1-based)}
//'   \item{N_samp_total}{\code{integer}: Total number of samples across all subjects}
//'   \item{samp_idx}{\code{integer} vector: Global sample index for each observation}
//'   \item{sample_to_subj_map}{\code{integer} vector: Maps global sample indices to subject indices}
//' }
//'
//' @details
//' The function performs several key transformations:
//' \itemize{
//'   \item Converts \code{character} subject IDs to consecutive numeric indices (1 to N_subj)
//'   \item Calculates the maximum number of samples per subject
//'   \item Creates global sample indexing that accounts for varying samples per subject
//'   \item Generates a mapping from global sample indices back to subject indices
//' }
//'
//' All input vectors must have identical lengths. The function includes comprehensive
//' error checking for empty data, mismatched vector lengths, and invalid sample IDs.
//' Sample IDs must be positive \code{integer}s starting from 1.
//'
//' @examples
//' \dontrun{
//' # Example with 3 subjects and varying samples per subject
//' subjects <- c("A", "A", "B", "B", "B", "C")
//' samples <- c(1, 2, 1, 2, 3, 1)
//' responses <- c(1.2, 1.5, 2.1, 2.3, 2.0, 0.8)
//' 
//' result <- process_stan_data_indexing(subjects, samples, responses)
//' print(result)
//' }
//'
//' @seealso
//' This function is designed for use with Stan models that require hierarchical
//' indexing structures for subjects and samples.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List process_stan_data_indexing(Rcpp::CharacterVector SubjectID_orig_R, 
                                      Rcpp::IntegerVector SampleID_orig_R, 
                                      Rcpp::NumericVector y_R) {
  int N_obs = SubjectID_orig_R.length();
  if (N_obs == 0) {
    Rcpp::stop("Input data is empty.");
  }
  if (SampleID_orig_R.length() != N_obs || y_R.length() != N_obs) {
    Rcpp::stop("Input vectors must have the same length.");
  }
  
  // Output vectors
  Rcpp::IntegerVector subj_idx(N_obs);
  Rcpp::IntegerVector samp_idx(N_obs);
  
  // 1. Subject mapping and N_subj
  std::map<std::string, int> subject_to_numeric_id;
  int N_subj = 0;
  
  for (int i = 0; i < N_obs; ++i) {
    std::string current_subj_str(SubjectID_orig_R[i]);
    if (subject_to_numeric_id.find(current_subj_str) == subject_to_numeric_id.end()) {
      N_subj++;
      subject_to_numeric_id[current_subj_str] = N_subj;
    }
    subj_idx[i] = subject_to_numeric_id[current_subj_str];
  }
  
  if (N_subj == 0 && N_obs > 0) {
    Rcpp::stop("N_subj is 0 but N_obs > 0. This should not happen.");
  }
  if (N_subj == 0 && N_obs == 0) { // Handle truly empty case if necessary
    return Rcpp::List::create(
      Rcpp::Named("N_obs") = 0,
      Rcpp::Named("y") = Rcpp::NumericVector(0),
      Rcpp::Named("N_subj") = 0,
      Rcpp::Named("subj_idx") = Rcpp::IntegerVector(0),
      Rcpp::Named("N_samp_total") = 0,
      Rcpp::Named("samp_idx") = Rcpp::IntegerVector(0),
      Rcpp::Named("sample_to_subj_map") = Rcpp::IntegerVector(0)
    );
  }
  
  
  // 2. Calculate max samples per subject (using 0-based numeric subject IDs internally for vector indexing)
  // max_samples_per_subj_numeric will be 0-indexed (for subject 1, index 0, etc.)
  std::vector<int> max_samples_per_subj_numeric(N_subj, 0); 
  for (int i = 0; i < N_obs; ++i) {
    int numeric_s_id = subj_idx[i] - 1; // 0-indexed for vector
    if (numeric_s_id < 0 || numeric_s_id >= N_subj) {
      Rcpp::Rcout << "Observation index: " << i << std::endl;
      Rcpp::Rcout << "subj_idx[i]: " << subj_idx[i] << std::endl;
      Rcpp::Rcout << "N_subj: " << N_subj << std::endl;
      Rcpp::stop("Error: numeric_s_id out of bounds during max_samples_per_subj_numeric calculation.");
    }
    max_samples_per_subj_numeric[numeric_s_id] = 
      std::max(max_samples_per_subj_numeric[numeric_s_id], SampleID_orig_R[i]);
  }
  
  // 3. N_samp_total: Total number of unique samples across all subjects
  int N_samp_total = 0;
  for (int count : max_samples_per_subj_numeric) {
    N_samp_total += count;
  }
  
  if (N_samp_total == 0 && N_obs > 0) {
    Rcpp::stop("N_samp_total is 0 but N_obs > 0. Check SampleID values (must be >0).");
  }
  
  
  // 4. Calculate cumulative offsets for global sample indexing (0-indexed offsets for 0-indexed subjects)
  std::vector<int> offsets_vec(N_subj);
  if (N_subj > 0) {
    offsets_vec[0] = 0;
    for (int s = 1; s < N_subj; ++s) {
      offsets_vec[s] = offsets_vec[s-1] + max_samples_per_subj_numeric[s-1];
    }
  }
  
  
  // 5. Calculate global sample index (samp_idx) for each observation
  for (int i = 0; i < N_obs; ++i) {
    int numeric_s_id = subj_idx[i] - 1; // 0-indexed
    if (numeric_s_id < 0 || numeric_s_id >= N_subj) {
      Rcpp::stop("Error: numeric_s_id out of bounds during samp_idx calculation.");
    }
    samp_idx[i] = offsets_vec[numeric_s_id] + SampleID_orig_R[i]; // SampleID_orig_R is 1-based
  }
  
  // 6. Create sample_to_subj_map: maps each global sample index (1..N_samp_total) to its subject index (1..N_subj)
  Rcpp::IntegerVector sample_to_subj_map(N_samp_total);
  int current_global_samp_idx = 0; // 0-indexed for vector access
  for (int s_numeric = 0; s_numeric < N_subj; ++s_numeric) { // s_numeric is 0-indexed subject ID
    int num_samples_for_this_subj = max_samples_per_subj_numeric[s_numeric];
    for (int local_samp_num = 1; local_samp_num <= num_samples_for_this_subj; ++local_samp_num) {
      if (current_global_samp_idx >= N_samp_total) {
        Rcpp::stop("Error constructing sample_to_subj_map: index out of bounds.");
      }
      sample_to_subj_map[current_global_samp_idx] = s_numeric + 1; // Store 1-based subject ID
      current_global_samp_idx++;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("N_obs") = N_obs,
    Rcpp::Named("y") = y_R, // Pass through y
    Rcpp::Named("N_subj") = N_subj,
    Rcpp::Named("subj_idx") = subj_idx,
    Rcpp::Named("N_samp_total") = N_samp_total,
    Rcpp::Named("samp_idx") = samp_idx,
    Rcpp::Named("sample_to_subj_map") = sample_to_subj_map
  );
}