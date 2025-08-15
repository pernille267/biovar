#include <Rcpp.h>
#include <numeric> // For std::iota
#include <algorithm> // For std::sort
#include <vector>
#include <cmath> // For std::pow, std::sqrt, ISNAN
using namespace Rcpp;

// Helper function for mean, skipping NAs
inline double mean_no_na(const Rcpp::NumericVector& x) {
  double sum = 0;
  int count = 0;
  for (int i = 0; i < x.size(); ++i) {
    if (!Rcpp::NumericVector::is_na(x[i])) {
      sum += x[i];
      count++;
    }
  }
  if (count == 0) return NA_REAL;
  return sum / count;
}

// Helper function for standard deviation, skipping NAs
inline double sd_no_na(const Rcpp::NumericVector& x) {
  double m = mean_no_na(x);
  if (ISNAN(m)) return NA_REAL;
  double sum_sq_diff = 0;
  int count = 0;
  for (int i = 0; i < x.size(); ++i) {
    if (!Rcpp::NumericVector::is_na(x[i])) {
      sum_sq_diff += std::pow(x[i] - m, 2.0);
      count++;
    }
  }
  if (count < 2) return NA_REAL; 
  return std::sqrt(sum_sq_diff / (count - 1));
}

// Subset samples for a particular subject
IntegerVector subset_samples(IntegerVector samples, IntegerVector subjects, int subject_id) {
  IntegerVector output;
  for (int i = 0; i < subjects.size(); i++) {
    if (subjects[i] == subject_id) {
      output.push_back(samples[i]);
    }
  }
  return output;
}

// Subset sample values for a particular subject
NumericVector subset_samples_values(NumericVector values, IntegerVector subjects, int subject_id) {
  NumericVector output;
  for (int i = 0; i < subjects.size(); i++) {
    if (subjects[i] == subject_id) {
      output.push_back(values[i]);
    }
  }
  return output;
}

// Subset replicates for a particular sample and subject
IntegerVector subset_replicates(IntegerVector replicates, IntegerVector subjects, IntegerVector samples, int subject_id, int sample_id) {
  IntegerVector output;
  for (int i = 0; i < replicates.size(); i++) {
    if (subjects[i] == subject_id && samples[i] == sample_id) {
      output.push_back(replicates[i]);
    }
  }
  return output;
}

// Subset replicate values for a particular sample and subject
NumericVector subset_replicates_values(NumericVector values, IntegerVector subjects, IntegerVector samples, int subject_id, int sample_id) {
  NumericVector output;
  for (int i = 0; i < values.size(); i++) {
    if (subjects[i] == subject_id && samples[i] == sample_id) {
      output.push_back(values[i]);
    }
  }
  return output;
}


//' One Way and Two Way ANOVA Modelling Based on Biological Variation Data
//'
//' @title One Way and Two Way ANOVA Modelling Based on Biological Variation Data
//' @name bv_anova
//'
//' @param data \code{data.table} or \code{list} object, with the following columns or elements:
//' \itemize{
//'   \item \code{SubjectID} - \code{Integer} vector with the IDs of the subjects.
//'   \item \code{SampleID} - \code{Integer} vector with the IDs of the samples.
//'   \item \code{ReplicateID} - \code{Integer} vector with the IDs of the replicates.
//'   \item \code{y} - \code{Numeric} vector with the measured values for each combination of IDs.
//' }
//' @param cv_anova \code{logical}. If \code{TRUE}, CV-ANOVA is performed.
//' 
//'
//' @description
//' This function uses biological variation data to estimate both one-way
//' (one model for each subject) and two-way (general model overall subjects)
//' nested ANOVA models. 
//' @details
//' The components of the estimated model can be used to estimate intra
//' biological variation, inter biological variation and analytical variation,
//' possibly with corresponding confidence intervals. The model estimates
//' relies on that all model effects are mutually independent and normally
//' distributed.
//' 
//' We assume the following model
//' \eqn{y_{isr} = \mu + G_i + I_{is} + A_{isr}}, where 
//' \eqn{G_i \sim \mathrm{N}(0, \sigma_G^2)},
//' \eqn{I_{is} \sim \mathrm{N}(0, \sigma_i^2)}, and
//' \eqn{A_{isr} \sim \mathrm{N}(0, \sigma_A^2)}. 
//'
//' @return A \code{list} containing all relevant information for the nested ANOVA model fit.
//'
//' @examples
//' \dontrun{
//' # Simulate data for 15 subjects, with 10 samples per subject and 3 replicates per sample
//' output <- simulate_bv_data(15, 10, 2, 10, 2, 50, 100)
//' anova_fit <- bv_anova(output)
//' print(anova_fit)
//' }
//'

// [[Rcpp::export]]
Rcpp::List bv_anova(Rcpp::List data, bool cv_anova = false) {
  
  Rcpp::IntegerVector subjects = data["SubjectID"];
  Rcpp::IntegerVector samples = data["SampleID"];
  Rcpp::NumericVector values_raw = data["y"];
  Rcpp::NumericVector values = clone(values_raw);
  
  int N_obs_total_rows = values.size(); // Total rows in input, might include NAs
  if (N_obs_total_rows == 0) {
    Rcpp::stop("Input data is empty.");
  }
  
  std::vector<int> p_idx(N_obs_total_rows);
  std::iota(p_idx.begin(), p_idx.end(), 0);
  
  std::sort(p_idx.begin(), p_idx.end(), [&](int a, int b) {
    if (subjects[a] != subjects[b]) {
      return subjects[a] < subjects[b];
    }
    if (samples[a] != samples[b]) {
      return samples[a] < samples[b];
    }
    return a < b; 
  });
  
  double grand_mean_val = mean_no_na(values);
  if (ISNAN(grand_mean_val)) {
    Rcpp::stop("Grand mean is missing or uncalculable (all NA values?).");
  }
  
  if (cv_anova) {
    values = values / grand_mean_val;
    grand_mean_val = 1.0;
  }
  
  double SST_val = 0;
  long long N_eff = 0; // Effective number of non-NA observations
  for (int k_obs = 0; k_obs < N_obs_total_rows; ++k_obs) {
    if (!Rcpp::NumericVector::is_na(values[k_obs])) {
      SST_val += std::pow(values[k_obs] - grand_mean_val, 2.0);
      N_eff++; // Count non-NA observations
    }
  }
  if (N_eff == 0) {
    Rcpp::stop("No non-NA observations found in 'y'.");
  }
  
  double SS1_val = 0, SS2_val = 0, SS3_val = 0;
  double SS1U_sum_sq_dev_grand_mean_acc = 0; // Sum of (subj_w_mean_i** - overall_w_grand_mean_***)^2
  double SS2U_sum_sq_dev_subj_mean_acc = 0;  // Sum over subjects of [ sum_j (sample_mean_ij. - subj_w_mean_i**)^2 ]
  
  int n_unique_subjects = 0;
  int n_total_unique_samples = 0;
  
  double w1U_num_acc = 0;       // For w1U: Sum_i (1/J_i)
  double w1U_denom_acc = 0;     // For w1U & w2U: Sum_i (X_i) = Sum_i (1 / (J_i * R_bar_i))
  double w2U_sum_X_i_sq_acc = 0; // For w2U: Sum_i (X_i^2)
  
  double w3U_textbook_num_acc = 0;  // For textbook w3U numerator: Sum_i ((J_i-1)/J_i)
  double w3U_textbook_den_acc = 0;  // For textbook w3U denominator: Sum_i ( ((J_i-1)/J_i) / R_bar_i )
  
  Rcpp::NumericVector SSi_per_subject_vec; 
  Rcpp::NumericVector SSa_per_subject_vec; 
  Rcpp::NumericVector Si_squared_per_subject_vec;
  Rcpp::NumericVector Sa_squared_per_subject_vec;
  Rcpp::IntegerVector ni_per_subject_vec; 
  Rcpp::IntegerVector na_per_subject_vec; 
  
  std::vector<double> subject_weighted_means_for_SS1U_vec; 
  double overall_weighted_grand_mean_numerator_acc = 0;
  
  int current_processing_idx = 0;
  while (current_processing_idx < N_obs_total_rows) {
    // Skip any leading NA values if the sorting placed them first
    while(current_processing_idx < N_obs_total_rows && Rcpp::NumericVector::is_na(values[p_idx[current_processing_idx]])) {
      current_processing_idx++;
    }
    if (current_processing_idx >= N_obs_total_rows) break;
    
    n_unique_subjects++;
    // int subject_loop_start_idx = current_processing_idx; // Not strictly used after NA skip
    int current_subject_original_id = subjects[p_idx[current_processing_idx]];
    
    double subject_total_y_sum = 0;    // Sum y_ijk for subject i
    int    subject_total_obs_count_Ni = 0; // N_i = Sum_j K_ij (non-NA obs for subject i)
    int    samples_in_this_subject_count_Ji = 0; // J_i (number of samples in subject i)
    double subject_sum_inv_Kij_for_Rbar = 0; // Sum_j (1/K_ij) for this subject
    double subject_sum_sample_means_for_ybar_i_starstar = 0; // Sum_j y_bar_ij. (for y_bar_i**)
    
    std::vector<double> y_values_this_subject_non_na; 
    std::vector<double> sample_means_this_subject_ybar_ij_dot; 
    std::vector<int>    replicates_per_sample_this_subject_Kij; 
    
    int subject_data_start_idx = current_processing_idx;
    
    while (current_processing_idx < N_obs_total_rows && subjects[p_idx[current_processing_idx]] == current_subject_original_id) {
      // Skip NA values within a subject's block before processing a sample
      while(current_processing_idx < N_obs_total_rows && 
            subjects[p_idx[current_processing_idx]] == current_subject_original_id &&
            Rcpp::NumericVector::is_na(values[p_idx[current_processing_idx]])) {
        current_processing_idx++;
      }
      if (current_processing_idx >= N_obs_total_rows || subjects[p_idx[current_processing_idx]] != current_subject_original_id) break;
      
      samples_in_this_subject_count_Ji++;
      n_total_unique_samples++;
      int sample_loop_start_idx = current_processing_idx;
      int current_sample_original_id = samples[p_idx[current_processing_idx]];
      
      double sample_y_sum = 0;
      int    sample_replicates_count_K_ij = 0;
      
      while (current_processing_idx < N_obs_total_rows &&
             subjects[p_idx[current_processing_idx]] == current_subject_original_id &&
             samples[p_idx[current_processing_idx]] == current_sample_original_id) {
        
        double y_val = values[p_idx[current_processing_idx]];
        if (!Rcpp::NumericVector::is_na(y_val)) {
          sample_y_sum += y_val;
          y_values_this_subject_non_na.push_back(y_val); 
          sample_replicates_count_K_ij++;
        }
        current_processing_idx++;
      }
      
      if (sample_replicates_count_K_ij > 0) {
        double sample_mean_ij_dot = sample_y_sum / sample_replicates_count_K_ij;
        sample_means_this_subject_ybar_ij_dot.push_back(sample_mean_ij_dot);
        replicates_per_sample_this_subject_Kij.push_back(sample_replicates_count_K_ij);
        
        subject_sum_sample_means_for_ybar_i_starstar += sample_mean_ij_dot;
        subject_sum_inv_Kij_for_Rbar += 1.0 / sample_replicates_count_K_ij;
        
        for (int k_repl_loop_idx = sample_loop_start_idx; k_repl_loop_idx < current_processing_idx; ++k_repl_loop_idx) {
          if (subjects[p_idx[k_repl_loop_idx]] == current_subject_original_id &&
              samples[p_idx[k_repl_loop_idx]] == current_sample_original_id && // Redundant if loop structure correct
              !Rcpp::NumericVector::is_na(values[p_idx[k_repl_loop_idx]])) {
              SS3_val += std::pow(values[p_idx[k_repl_loop_idx]] - sample_mean_ij_dot, 2.0);
          }
        }
        subject_total_y_sum += sample_y_sum;
        subject_total_obs_count_Ni += sample_replicates_count_K_ij;
      } else {
        // If a sample ID was present but all its values were NA, it doesn't contribute to K_ij counts
        // but might have incremented samples_in_this_subject_count_Ji if it was a distinct sample ID.
        // For Ji calculation, it should be samples with at least one non-NA observation.
        // Let's adjust Ji if sample_replicates_count_K_ij was 0.
        samples_in_this_subject_count_Ji--; 
        n_total_unique_samples--;
      }
    } // End of Samples for this Subject
    
    if (subject_total_obs_count_Ni > 0 && samples_in_this_subject_count_Ji > 0) {
      double subject_mean_i_dotdot = subject_total_y_sum / subject_total_obs_count_Ni; // y_bar_i..
      SS1_val += subject_total_obs_count_Ni * std::pow(subject_mean_i_dotdot - grand_mean_val, 2.0);
      
      double current_SSt_i = 0;
      for(double val_y_subj : y_values_this_subject_non_na) {
        current_SSt_i += std::pow(val_y_subj - subject_mean_i_dotdot, 2.0);
      }
      
      double current_SSi_i = 0; // SSi for subject i (SSB(i))
      for (size_t j = 0; j < sample_means_this_subject_ybar_ij_dot.size(); ++j) {
        current_SSi_i += replicates_per_sample_this_subject_Kij[j] * 
          std::pow(sample_means_this_subject_ybar_ij_dot[j] - subject_mean_i_dotdot, 2.0);
      }
      SS2_val += current_SSi_i;
      SSi_per_subject_vec.push_back(current_SSi_i);
      
      double current_SSa_i = current_SSt_i - current_SSi_i; // SSa for subject i (SSE(i))
      SSa_per_subject_vec.push_back(current_SSa_i);
      
      int current_ni_df_Ji_minus_1 = samples_in_this_subject_count_Ji - 1;
      int current_na_df_Ni_minus_Ji = subject_total_obs_count_Ni - samples_in_this_subject_count_Ji;
      ni_per_subject_vec.push_back(current_ni_df_Ji_minus_1);
      na_per_subject_vec.push_back(current_na_df_Ni_minus_Ji);
      
      Si_squared_per_subject_vec.push_back((current_ni_df_Ji_minus_1 > 0) ? current_SSi_i / current_ni_df_Ji_minus_1 : NA_REAL);
      Sa_squared_per_subject_vec.push_back((current_na_df_Ni_minus_Ji > 0) ? current_SSa_i / current_na_df_Ni_minus_Ji : NA_REAL);
      
      // For SS1U, SS2U, and w coefficients (using textbook definitions)
      // y_bar_i** (textbook notation) = (1/J_i) * Sum_j y_bar_ij.
      double subject_weighted_mean_ybar_i_starstar = subject_sum_sample_means_for_ybar_i_starstar / samples_in_this_subject_count_Ji;
      subject_weighted_means_for_SS1U_vec.push_back(subject_weighted_mean_ybar_i_starstar);
      overall_weighted_grand_mean_numerator_acc += subject_weighted_mean_ybar_i_starstar;
      
      double ss2u_contrib_this_subj = 0;
      for(double sm_ij_dot : sample_means_this_subject_ybar_ij_dot) {
        ss2u_contrib_this_subj += std::pow(sm_ij_dot - subject_weighted_mean_ybar_i_starstar, 2.0);
      }
      SS2U_sum_sq_dev_subj_mean_acc += ss2u_contrib_this_subj;
      
      // R_bar_i. = J_i / (Sum_j 1/K_ij)
      double R_bar_i_dot = 0.0;
      if (subject_sum_inv_Kij_for_Rbar > 0) {
        R_bar_i_dot = static_cast<double>(samples_in_this_subject_count_Ji) / subject_sum_inv_Kij_for_Rbar;
      }
      
      if (R_bar_i_dot > 0) { // Avoid division by zero if R_bar_i_dot is 0 (e.g. all K_ij were 0, though filtered)
        double J_i_val = static_cast<double>(samples_in_this_subject_count_Ji);
        
        w1U_num_acc += 1.0 / J_i_val;
        double X_i = 1.0 / (J_i_val * R_bar_i_dot);
        w1U_denom_acc += X_i; // This is Sum_i X_i
        w2U_sum_X_i_sq_acc += X_i * X_i; // Sum_i X_i^2
        
        if (J_i_val > 0) { // Should always be true here due to outer if
          double Y_i_term = (J_i_val - 1.0) / J_i_val;
          w3U_textbook_num_acc += Y_i_term;
          w3U_textbook_den_acc += Y_i_term / R_bar_i_dot;
        }
      }
    } else if (n_unique_subjects > 0 && subject_total_obs_count_Ni == 0 && samples_in_this_subject_count_Ji == 0) {
      // This case means a subject ID was encountered, but all its observations were NA.
      // Decrement n_unique_subjects because it's not contributing to ANOVA.
      n_unique_subjects--;
    }
  } // End of Main Subject Loop
  
  if (n_unique_subjects == 0) {
    Rcpp::warning("No subjects with valid data found after NA processing.");
    // Return a list of NAs or empty structures
    return Rcpp::List::create(Rcpp::Named("Error") = "No valid subjects after NA processing.");
  }
  
  int df1_n1 = (n_unique_subjects > 0) ? n_unique_subjects - 1 : 0;
  int df2_n2 = n_total_unique_samples - n_unique_subjects; 
  // df2_n2 must be non-negative. If n_total_unique_samples < n_unique_subjects (should not happen with correct logic), set to 0.
  if (df2_n2 < 0) df2_n2 = 0;
  
  int df3_n3 = (N_eff >= static_cast<long long>(n_total_unique_samples)) ? (static_cast<int>(N_eff) - n_total_unique_samples) : 0;
  if (df3_n3 < 0) df3_n3 = 0; // Should not happen if N_eff, n_total_unique_samples are correct
  
  double overall_weighted_grand_mean_ybar_starstarstar = 0;
  if (n_unique_subjects > 0) {
    overall_weighted_grand_mean_ybar_starstarstar = overall_weighted_grand_mean_numerator_acc / n_unique_subjects;
  }
  
  double final_SS1U_sum_sq_dev_num = 0; // Sum_i (y_bar_i** - y_bar_***)^2
  for(double swm_i_starstar : subject_weighted_means_for_SS1U_vec) {
    final_SS1U_sum_sq_dev_num += std::pow(swm_i_starstar - overall_weighted_grand_mean_ybar_starstarstar, 2.0);
  }
  
  double w1U = NA_REAL, w2U = NA_REAL, w3U = NA_REAL;
  double I_val = static_cast<double>(n_unique_subjects);
  
  if (w1U_denom_acc != 0.0) { // w1U_denom_acc is Sum_i X_i
    w1U = w1U_num_acc / w1U_denom_acc;
  }
  
  if (df1_n1 > 0 && I_val > 0 && w1U_denom_acc != 0.0) { // df1_n1 = I-1
    double sum_X_i = w1U_denom_acc;
    double sum_X_i_sq = w2U_sum_X_i_sq_acc;
    double term_in_bracket_num_for_w2U = I_val - (sum_X_i_sq / (sum_X_i * sum_X_i)); // Handles sum_X_i * sum_X_i = 0 if sum_X_i is 0
    double term_denom_for_w2U = (1.0 / I_val) * sum_X_i;
    if (term_denom_for_w2U != 0.0) {
      w2U = (1.0 / static_cast<double>(df1_n1)) * term_in_bracket_num_for_w2U / term_denom_for_w2U;
    }
  }
  
  if (df2_n2 > 0 && w3U_textbook_den_acc != 0.0) {
    w3U = w3U_textbook_num_acc / w3U_textbook_den_acc;
  } else if (df2_n2 > 0 && w3U_textbook_num_acc == 0.0 && w3U_textbook_den_acc == 0.0) {
    w3U = 0.0; // Case where all J_i=1 -> num=0, den=0, so SS2U will be 0. Or NA_REAL if preferred.
  }
  
  double final_SS1U = NA_REAL, final_SS2U = NA_REAL;
  if (!ISNAN(w2U)) final_SS1U = final_SS1U_sum_sq_dev_num * w2U;
  if (!ISNAN(w3U)) final_SS2U = SS2U_sum_sq_dev_subj_mean_acc * w3U;
  
  double S1_squared  = (df1_n1 > 0) ? SS1_val / df1_n1 : NA_REAL;
  double S2_squared  = (df2_n2 > 0) ? SS2_val / df2_n2 : NA_REAL;
  double S3_squared  = (df3_n3 > 0) ? SS3_val / df3_n3 : NA_REAL;
  
  double S1U_squared = (df1_n1 > 0 && !ISNAN(final_SS1U)) ? final_SS1U / df1_n1 : NA_REAL;
  double S2U_squared = (df2_n2 > 0 && !ISNAN(final_SS2U)) ? final_SS2U / df2_n2 : NA_REAL;
  
  Rcpp::NumericVector Sigma_i_vec(n_unique_subjects);
  Rcpp::NumericVector Sigma_a_vec(n_unique_subjects);
  
  for(int i=0; i < n_unique_subjects; ++i) { // n_unique_subjects is now the count of subjects with data
    if (i < Sa_squared_per_subject_vec.size() && !ISNAN(Sa_squared_per_subject_vec[i])) {
      Sigma_a_vec[i] = std::sqrt(Sa_squared_per_subject_vec[i]);
    } else {
      Sigma_a_vec[i] = NA_REAL;
    }
    
    if (i < Si_squared_per_subject_vec.size() && i < Sa_squared_per_subject_vec.size() &&
        !ISNAN(Si_squared_per_subject_vec[i]) && !ISNAN(Sa_squared_per_subject_vec[i]) && 
        !ISNAN(w3U) && w3U != 0.0) {
        double term_under_sqrt = (Si_squared_per_subject_vec[i] - Sa_squared_per_subject_vec[i]) / w3U;
      // Original formulation was Si_squared/w3U - Sa_squared/w3U.
      // If EMS_B(i) = sigma_e_sq(i) + k_i * sigma_b_sq(i), then sigma_b_sq(i) = (MSB(i) - MSE(i)) / k_i
      // Here it seems w3U is used as a global k_i proxy.
      Sigma_i_vec[i] = (term_under_sqrt >= 0) ? std::sqrt(term_under_sqrt) : NA_REAL; 
    } else {
      Sigma_i_vec[i] = NA_REAL;
    }
  }
  
  double hbhr_val = NA_REAL;
  if (Sigma_i_vec.size() > 0) { // Ensure Sigma_i_vec is not empty
    double mean_sigma_i = mean_no_na(Sigma_i_vec);
    if (!ISNAN(mean_sigma_i) && mean_sigma_i != 0) { // Avoid division by zero
      double sd_sigma_i = sd_no_na(Sigma_i_vec);
      if(!ISNAN(sd_sigma_i)) {
        hbhr_val = sd_sigma_i / mean_sigma_i * 100.0;
      }
    } else if (!ISNAN(mean_sigma_i) && mean_sigma_i == 0) {
      // If mean is 0, check if sd is also 0
      double sd_sigma_i = sd_no_na(Sigma_i_vec);
      if (!ISNAN(sd_sigma_i) && sd_sigma_i == 0.0) {
        hbhr_val = 0.0; // All sigma_i are 0
      }
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("grand_mean") = grand_mean_val,
    Rcpp::Named("S1_squared") = S1_squared,
    Rcpp::Named("S2_squared") = S2_squared,
    Rcpp::Named("S3_squared") = S3_squared,
    Rcpp::Named("S1U_squared") = S1U_squared,
    Rcpp::Named("S2U_squared") = S2U_squared,
    Rcpp::Named("n1") = df1_n1, 
    Rcpp::Named("n2") = df2_n2,
    Rcpp::Named("n3") = df3_n3,
    Rcpp::Named("N_eff_obs") = N_eff,
    Rcpp::Named("I_subjects") = n_unique_subjects,
    Rcpp::Named("J_total_samples") = n_total_unique_samples,
    Rcpp::Named("w1U") = w1U,
    Rcpp::Named("w2U") = w2U,
    Rcpp::Named("w3U") = w3U,
    Rcpp::Named("SS1U") = final_SS1U,
    Rcpp::Named("SS2U") = final_SS2U,
    Rcpp::Named("Si_squared_persubject") = Si_squared_per_subject_vec, 
    Rcpp::Named("Sa_squared_persubject") = Sa_squared_per_subject_vec, 
    Rcpp::Named("Sigma_i") = Sigma_i_vec,                   
    Rcpp::Named("Sigma_a") = Sigma_a_vec,                   
    Rcpp::Named("ni") = ni_per_subject_vec,                 
    Rcpp::Named("na") = na_per_subject_vec,                 
    Rcpp::Named("hbhr") = hbhr_val,
    Rcpp::Named("SST") = SST_val,
    Rcpp::Named("SS1_FactorA") = SS1_val,
    Rcpp::Named("SS2_FactorB_within_A") = SS2_val,
    Rcpp::Named("SS3_Error") = SS3_val
  );
}

//' @title Estimate Variance Components and Confidence Intervals via ANOVA
//'
//' @name variance_components
//'
//' @param data A \code{list} or \code{data.frame} containing the data for the ANOVA.
//'   It must contain columns for subject identifiers and measurement values. This
//'   is passed to an internal ANOVA function (`bv_anova`).
//' @param output_type A \code{character} string specifying the type of output to return.
//'   Must be one of the following:
//'   \describe{
//'     \item{\code{sigma}}{Returns point estimates of the standard deviations (\code{sigma}).}
//'     \item{\code{cv}}{Returns point estimates of the coefficients of variation (CV), calculated as sigma / grand_mean.}
//'     \item{\code{sigma_ci}}{Returns point estimates and confidence intervals for the standard deviations.}
//'     \item{\code{cv_ci}}{Returns point estimates and confidence intervals for the CVs.}
//'   }
//'   Defaults to \code{"sigma"}.
//' @param mult A \code{double} used as a scaling factor for the final output. For
//'   example, set to 100 to express CVs as percentages. Defaults to 1.0.
//' @param level A \code{double} specifying the confidence level for the confidence
//'   intervals (e.g., 0.95 for a 95\% CI). Defaults to 0.95.
//' @param cv_anova A \code{bool}. If \code{TRUE}, the underlying ANOVA is performed on
//'   subject-specific coefficients of variation rather than on the raw
//'   measurement values. Defaults to \code{FALSE}.
//'
//' @description
//' This function performs a three-level nested ANOVA to estimate variance
//' components. It can return point estimates of standard deviations (sigmas),
//' coefficients of variation (CVs), or confidence intervals for these metrics.
//'
//' @details
//' The function is built on a three-level nested random-effects model, which
//' decomposes the total variance into three components:
//' \itemize{
//'   \item \strong{\code{sigma_G}}: The top-level, between-subject standard deviation.
//'   \item \strong{\code{sigma_I}}: The intermediate-level, within-subject standard deviation.
//'   \item \strong{\code{sigma_A}}: The lowest-level, within-replicate (analytical or residual) standard deviation.
//' }
//' Point estimates for these components are derived from the mean squares of the
//' ANOVA. Confidence intervals are also computed. The CI for the residual variance
//' (\code{sigma_A}) is based on the chi-squared distribution. For the higher-level
//' components (\code{sigma_I} and \code{sigma_G}), more complex approximate methods
//' (generalized confidence intervals) are used, which are suitable for unbalanced data.
//' The function also returns subject-specific point estimates and the
//' Heterogeneity of the Biological Homeostatic Ratio (HBHR).
//'
//' @return A \code{list} containing the estimated components. The structure of
//'   the first three elements depends on \code{output_type}.
//'   \describe{
//'     \item{\code{sigma_A}}{Analytical/residual component. A single \code{double} for point
//'       estimates, or a 3-element \code{NumericVector} (estimate, lower CI, upper CI) for CIs.}
//'     \item{\code{sigma_I}}{Within-subject component. A single \code{double} or a 3-element \code{NumericVector}.}
//'     \item{\code{sigma_G}}{Between-subject component. A single \code{double} or a 3-element \code{NumericVector}.}
//'     \item{\code{sigma_i}}{A \code{NumericVector} of point estimates for the subject-specific within-subject standard deviation.}
//'     \item{\code{sigma_a}}{A \code{NumericVector} of point estimates for the subject-specific analytical standard deviation.}
//'     \item{\code{HBHR}}{A \code{double} representing the Heterogeneity of the Biological Homeostatic Ratio.}
//'   }
//' @export
//' @examples
//' # Create a sample data frame for demonstration
//' # 3 subjects, 3 time points per subject, 2 replicates per time point
//' subject_data <- data.frame(
//'   SubjectID = rep(c("S1", "S2", "S3"), each = 6),
//'   Value = c(
//'     rnorm(2, 100, 10), rnorm(2, 105, 10), rnorm(2, 98, 10),  # Subject 1
//'     rnorm(2, 120, 12), rnorm(2, 118, 12), rnorm(2, 122, 12), # Subject 2
//'     rnorm(2, 90, 8), rnorm(2, 92, 8), rnorm(2, 88, 8)       # Subject 3
//'   )
//' )
//'
//' # Note: In a real package, the bv_anova function would need to be available.
//' # The following calls are illustrative of how variance_components would be used.
//'
//' # Get point estimates of standard deviations (sigmas)
//' # results_sigma <- variance_components(subject_data, output_type = "sigma")
//'
//' # Get CVs as percentages with 95% confidence intervals
//' # results_cv_ci <- variance_components(
//' #   subject_data,
//' #   output_type = "cv_ci",
//' #   mult = 100,
//' #   level = 0.95
//' # )
//'
//' # print(results_sigma)
//' # print(results_cv_ci)
// [[Rcpp::export]]
Rcpp::List variance_components(List data, 
                               std::string output_type = "sigma", 
                               double mult = 1.0, 
                               double level = 0.95,
                               bool cv_anova = false) {
  
  // Run ANOVA model
  List anova_output_list = bv_anova(data, cv_anova);
  
  // Extract values from the computed ANOVA list
  double S1_squared = Rcpp::as<double>(anova_output_list["S1U_squared"]);
  double S2_squared = Rcpp::as<double>(anova_output_list["S2U_squared"]);
  double S3_squared = Rcpp::as<double>(anova_output_list["S3_squared"]);
  double w1U = Rcpp::as<double>(anova_output_list["w1U"]);
  double w2U = Rcpp::as<double>(anova_output_list["w2U"]);
  double w3U = Rcpp::as<double>(anova_output_list["w3U"]);
  double grand_mean = Rcpp::as<double>(anova_output_list["grand_mean"]);
  Rcpp::NumericVector sigma_i_persubj_point = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(anova_output_list["Sigma_i"]));
  Rcpp::NumericVector sigma_a_persubj_point = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(anova_output_list["Sigma_a"]));
  double HBHR = Rcpp::as<double>(anova_output_list["hbhr"]);
  
  // Point estimates for variance components (sigmas, not variances initially)
  double sigma_A_point = NA_REAL; // sigma_E (Analytical error)
  double sigma_I_point = NA_REAL; // sigma_B (within-subject)
  double sigma_G_point = NA_REAL; // sigma_A (between-subjects)
  
  // Initialize ratios c3 = w1u / w3u and c2 = c3 - 1. c3 = 1 if design is balanced.
  double c2 = NA_REAL, c3 = NA_REAL;
  
  // Calculate sigma_A_point (sigma_E)
  if (!ISNAN(S3_squared) && S3_squared >= 0) {
    sigma_A_point = std::sqrt(S3_squared);
  }
  
  // Calculate c2, c3, and sigma_I_point (sigma_B)
  if (!ISNAN(w1U) && !ISNAN(w3U) && w3U != 0) {
    c2 = w1U / w3U;
    c3 = c2 - 1.0; // c3 can be negative
    if (!ISNAN(S2_squared) && !ISNAN(S3_squared)) {
      double var_I_hat = (S2_squared - S3_squared) / w3U;
      if (var_I_hat >= 0) {
        sigma_I_point = std::sqrt(var_I_hat);
      } else {
        sigma_I_point = 0;
      }
    }
  }
  
  // Calculate sigma_G_point (sigma_A - top level)
  if (!ISNAN(w2U) && w2U != 0 && !ISNAN(c2) && !ISNAN(c3) &&
      !ISNAN(S1_squared) && !ISNAN(S2_squared) && !ISNAN(S3_squared)) {
      double var_G_hat = (S1_squared - c2 * S2_squared + c3 * S3_squared) / w2U;
    if (var_G_hat >= 0) {
      sigma_G_point = std::sqrt(var_G_hat);
    } else {
      sigma_G_point = 0; // Or NA_REAL
    }
  }
  
  // --- Handle "sigma" and "cv" outputs (without CIs) ---
  if (output_type == "sigma") {
    for (int i = 0; i < sigma_i_persubj_point.size(); ++i) {
      if(!ISNAN(sigma_i_persubj_point[i])) sigma_i_persubj_point[i] *= mult;
      if(!ISNAN(sigma_a_persubj_point[i])) sigma_a_persubj_point[i] *= mult;
    }
    return Rcpp::List::create(Rcpp::Named("sigma_A") = ISNAN(sigma_A_point) ? NA_REAL : sigma_A_point * mult,
                              Rcpp::Named("sigma_I") = ISNAN(sigma_I_point) ? NA_REAL : sigma_I_point * mult,
                              Rcpp::Named("sigma_G") = ISNAN(sigma_G_point) ? NA_REAL : sigma_G_point * mult,
                              Rcpp::Named("sigma_i") = sigma_i_persubj_point,
                              Rcpp::Named("sigma_a") = sigma_a_persubj_point,
                              Rcpp::Named("HBHR") = HBHR,
                              Rcpp::Named("beta") = grand_mean);  
  } else if (output_type == "cv") {
    if (ISNAN(grand_mean) || grand_mean == 0) {
      Rcpp::warning("Grand mean is NA or zero; CVs cannot be calculated and will be NA.");
      for (int i = 0; i < sigma_i_persubj_point.size(); ++i) {
        sigma_i_persubj_point[i] = NA_REAL;
        sigma_a_persubj_point[i] = NA_REAL;
      }
      return Rcpp::List::create(Rcpp::Named("sigma_A") = NA_REAL,
                                Rcpp::Named("sigma_I") = NA_REAL,
                                Rcpp::Named("sigma_G") = NA_REAL,
                                Rcpp::Named("sigma_i") = sigma_i_persubj_point,
                                Rcpp::Named("sigma_a") = sigma_a_persubj_point,
                                Rcpp::Named("HBHR") = HBHR,
                                Rcpp::Named("beta") = grand_mean);
    }
    for (int i = 0; i < sigma_i_persubj_point.size(); ++i) {
      if(!ISNAN(sigma_i_persubj_point[i])) sigma_i_persubj_point[i] = sigma_i_persubj_point[i] / grand_mean * mult;
      if(!ISNAN(sigma_a_persubj_point[i])) sigma_a_persubj_point[i] = sigma_a_persubj_point[i] / grand_mean * mult;
    }
    return Rcpp::List::create(Rcpp::Named("sigma_A") = ISNAN(sigma_A_point) ? NA_REAL : sigma_A_point / grand_mean * mult,
                              Rcpp::Named("sigma_I") = ISNAN(sigma_I_point) ? NA_REAL : sigma_I_point / grand_mean * mult,
                              Rcpp::Named("sigma_G") = ISNAN(sigma_G_point) ? NA_REAL : sigma_G_point / grand_mean * mult,
                              Rcpp::Named("sigma_i") = sigma_i_persubj_point,
                              Rcpp::Named("sigma_a") = sigma_a_persubj_point,
                              Rcpp::Named("HBHR") = HBHR,
                              Rcpp::Named("beta") = grand_mean);
  }
  
  // --- Calculate Confidence Intervals (only if output_type is "sigma_ci" or "cv_ci") ---
  if (output_type != "sigma_ci" && output_type != "cv_ci") {
    Rcpp::stop("Invalid output_type. Must be 'sigma', 'cv', 'sigma_ci', or 'cv_ci'.");
  }
  
  // Extract degrees of freedom
  double n1 = Rcpp::as<double>(anova_output_list["n1"]); // Degrees of freedom
  double n2 = Rcpp::as<double>(anova_output_list["n2"]);
  double n3 = Rcpp::as<double>(anova_output_list["n3"]);
  
  double alpha_param = (1.0 - level) / 2.0; // For two-sided CIs
  
  Rcpp::NumericVector sigma_A_ci_vec(3); sigma_A_ci_vec[0] = sigma_A_point; sigma_A_ci_vec[1]=NA_REAL; sigma_A_ci_vec[2]=NA_REAL;
  Rcpp::NumericVector sigma_I_ci_vec(3); sigma_I_ci_vec[0] = sigma_I_point; sigma_I_ci_vec[1]=NA_REAL; sigma_I_ci_vec[2]=NA_REAL;
  Rcpp::NumericVector sigma_G_ci_vec(3); sigma_G_ci_vec[0] = sigma_G_point; sigma_G_ci_vec[1]=NA_REAL; sigma_G_ci_vec[2]=NA_REAL;
  
  // CI for sigma_A (Error component sigma_e)
  if (n3 > 0 && !ISNAN(S3_squared) && S3_squared >= 0) {
    double chi_sq_lower_quantile = R::qchisq(alpha_param, n3, true, false); // for upper bound of variance
    double chi_sq_upper_quantile = R::qchisq(1.0 - alpha_param, n3, true, false); // for lower bound of variance
    
    if (chi_sq_lower_quantile > 0 && !ISNAN(chi_sq_lower_quantile)) {
      sigma_A_ci_vec[2] = std::sqrt((n3 * S3_squared) / chi_sq_lower_quantile); // Upper CI limit
    }
    if (chi_sq_upper_quantile > 0 && !ISNAN(chi_sq_upper_quantile)) {
      sigma_A_ci_vec[1] = std::sqrt((n3 * S3_squared) / chi_sq_upper_quantile); // Lower CI limit
    }
    if (!ISNAN(sigma_A_ci_vec[1]) && !ISNAN(sigma_A_ci_vec[2]) && sigma_A_ci_vec[1] > sigma_A_ci_vec[2]) {
      std::swap(sigma_A_ci_vec[1], sigma_A_ci_vec[2]); // Ensure lwr < upr
    }
    if (ISNAN(sigma_A_ci_vec[1])) sigma_A_ci_vec[1] = 0; // Truncate at 0 if lower bound is problematic
  }
  
  // G and H coefficients
  double G1=NA_REAL, G2=NA_REAL, G3=NA_REAL, H1=NA_REAL, H2=NA_REAL, H3=NA_REAL;
  double G12=NA_REAL, G13=NA_REAL, G13_star=NA_REAL, G23=NA_REAL, G32=NA_REAL;
  double H12=NA_REAL, H13=NA_REAL, H23_star=NA_REAL, H23=NA_REAL, H32=NA_REAL; // Renamed H23 to H23
  
  bool can_calc_coeffs = (n1 > 0 && n2 > 0 && n3 > 0 && 
                          !ISNAN(S1_squared) && !ISNAN(S2_squared) && !ISNAN(S3_squared) &&
                          S1_squared >=0 && S2_squared >=0 && S3_squared >=0);
  
  if (can_calc_coeffs) {
    double F_a_n1_inf = R::qf(alpha_param, n1, 1e12, true, false);
    double F_a_n2_inf = R::qf(alpha_param, n2, 1e12, true, false);
    double F_a_n3_inf = R::qf(alpha_param, n3, 1e12, true, false);
    double F_1a_n1_inf = R::qf(1.0 - alpha_param, n1, 1e12, true, false);
    double F_1a_n2_inf = R::qf(1.0 - alpha_param, n2, 1e12, true, false);
    double F_1a_n3_inf = R::qf(1.0 - alpha_param, n3, 1e12, true, false);
    
    if (F_a_n1_inf > 0 && !ISNAN(F_a_n1_inf)) G1 = 1.0 - 1.0 / F_a_n1_inf; else G1 = 0; 
    if (F_a_n2_inf > 0 && !ISNAN(F_a_n2_inf)) G2 = 1.0 - 1.0 / F_a_n2_inf; else G2 = 0;
    if (F_a_n3_inf > 0 && !ISNAN(F_a_n3_inf)) G3 = 1.0 - 1.0 / F_a_n3_inf; else G3 = 0;
    if (F_1a_n1_inf > 0 && !ISNAN(F_1a_n1_inf)) H1 = 1.0 / F_1a_n1_inf - 1.0; else H1 = 0;
    if (F_1a_n2_inf > 0 && !ISNAN(F_1a_n2_inf)) H2 = 1.0 / F_1a_n2_inf - 1.0; else H2 = 0;
    if (F_1a_n3_inf > 0 && !ISNAN(F_1a_n3_inf)) H3 = 1.0 / F_1a_n3_inf - 1.0; else H3 = 0;
    
    double F_a_n1_n2 = R::qf(alpha_param, n1, n2, true, false);
    double F_a_n1_n3 = R::qf(alpha_param, n1, n3, true, false);
    double F_a_n2_n3 = R::qf(alpha_param, n2, n3, true, false);
    double F_a_n3_n2 = R::qf(alpha_param, n3, n2, true, false);
    double F_1a_n1_n2 = R::qf(1.0 - alpha_param, n1, n2, true, false);
    double F_1a_n1_n3 = R::qf(1.0 - alpha_param, n1, n3, true, false);
    double F_1a_n2_n3 = R::qf(1.0 - alpha_param, n2, n3, true, false);
    double F_1a_n3_n2 = R::qf(1.0 - alpha_param, n3, n2, true, false);
    double F_a_n1n3_inf = R::qf(alpha_param, n1 + n3, 1e12, true, false);
    double F_a_n2n3_inf = R::qf(alpha_param, n2 + n3, 1e12, true, false);
    
    // G_12, G_13 and H_12, H_13 
    if (!ISNAN(G1)&&!ISNAN(H2)&&!ISNAN(F_a_n1_n2)&&F_a_n1_n2>0) G12 = (std::pow(F_a_n1_n2 - 1.0, 2.0) - std::pow(G1, 2.0) * std::pow(F_a_n1_n2, 2.0) - std::pow(H2,2.0)) / F_a_n1_n2;
    if (!ISNAN(G1)&&!ISNAN(H3)&&!ISNAN(F_a_n1_n3)&&F_a_n1_n3>0) G13 = (std::pow(F_a_n1_n3 - 1.0, 2.0) - std::pow(G1, 2.0) * std::pow(F_a_n1_n3, 2.0) - std::pow(H3,2.0)) / F_a_n1_n3;
    if (!ISNAN(H1)&&!ISNAN(G2)&&!ISNAN(F_1a_n1_n2)&&F_1a_n1_n2>0) H12 = (std::pow(1.0 - F_1a_n1_n2, 2.0) - std::pow(H1, 2.0) * std::pow(F_1a_n1_n2, 2.0) - std::pow(G2, 2.0)) /F_1a_n1_n2;
    if (!ISNAN(H1)&&!ISNAN(G3)&&!ISNAN(F_1a_n1_n3)&&F_1a_n1_n3>0) H13 = (std::pow(1.0 - F_1a_n1_n3, 2.0) - std::pow(H1, 2.0) * std::pow(F_1a_n1_n3, 2.0) - std::pow(G3, 2.0)) /F_1a_n1_n3;
    // G_23, H_23, G32, H32
    if (!ISNAN(G2)&&!ISNAN(H3)&&!ISNAN(F_a_n2_n3)&&F_a_n2_n3>0) G23 = (std::pow(F_a_n2_n3 - 1.0, 2.0) - std::pow(G2, 2.0) * std::pow(F_a_n2_n3, 2.0) - std::pow(H3, 2.0)) / F_a_n2_n3;
    if (!ISNAN(H2)&&!ISNAN(G3)&&!ISNAN(F_1a_n2_n3)&&F_1a_n2_n3>0) H23 = (std::pow(1.0 - F_1a_n2_n3, 2.0) - std::pow(H2, 2.0) * std::pow(F_1a_n2_n3, 2.0) - std::pow(G3, 2.0)) / F_1a_n2_n3;
    if (!ISNAN(G3)&&!ISNAN(H2)&&!ISNAN(F_a_n3_n2)&&F_a_n3_n2>0) G32 = (std::pow(F_a_n3_n2 - 1.0, 2.0) - std::pow(G3, 2.0) * std::pow(F_a_n3_n2, 2.0) - std::pow(H2, 2.0)) / F_a_n3_n2;
    if (!ISNAN(H3)&&!ISNAN(G2)&&!ISNAN(F_1a_n3_n2)&&F_1a_n3_n2>0) H32 = (std::pow(1.0 - F_1a_n3_n2, 2.0) - std::pow(H3, 2.0) * std::pow(F_1a_n3_n2, 2.0) - std::pow(G2, 2.0)) / F_1a_n3_n2;
    // G_13^star, H_13^star
    if (!ISNAN(G1)&&!ISNAN(G3)&&!ISNAN(F_a_n1n3_inf)&&F_a_n1n3_inf>0&&n1>0&&n3>0)
      G13_star = std::pow(1.0 - 1.0 / F_a_n1n3_inf, 2.0) * std::pow(n1 + n3, 2.0) / (n1 * n3) - std::pow(G1, 2.0) * n1 / n3 - std::pow(G3, 2.0)* n3 / n1;
    if (!ISNAN(G2)&&!ISNAN(G3)&&!ISNAN(F_a_n2n3_inf)&&F_a_n2n3_inf>0&&n2>0&&n3>0)
      H23_star = std::pow(1.0 - 1.0 / F_a_n2n3_inf, 2.0) * std::pow(n2 + n3, 2.0) / (n2 * n3) - std::pow(G2, 2.0)* n2 / n3 - std::pow(G3, 2.0)* n3 / n2;
    
    // Confidence interval for sigma_I (sigma_b)
    if (!ISNAN(w3U) && w3U != 0 &&
        !ISNAN(G2) && !ISNAN(H3) && !ISNAN(G23) && 
        !ISNAN(H2) && !ISNAN(G3) && !ISNAN(H23)) {
      double V_L_I = std::pow(G2, 2.0) * std::pow(S2_squared, 2.0) + std::pow(H3, 2.0) * std::pow(S3_squared, 2.0) + G23 * S2_squared * S3_squared;
      double V_U_I = std::pow(H2, 2.0) * std::pow(S2_squared, 2.0) + std::pow(G3, 2.0) * std::pow(S3_squared, 2.0) + H23 * S2_squared * S3_squared;
      if (V_L_I >= 0 && V_U_I >= 0) {
        double term_I_var_hat = (S2_squared - S3_squared) / w3U;
        double lower_var_I = term_I_var_hat - std::sqrt(V_L_I) / w3U; 
        double upper_var_I = term_I_var_hat + std::sqrt(V_U_I) / w3U;
        
        sigma_I_ci_vec[1] = std::sqrt(std::max(0.0, lower_var_I));
        sigma_I_ci_vec[2] = std::sqrt(std::max(0.0, upper_var_I));
        if (!ISNAN(sigma_I_ci_vec[1]) && !ISNAN(sigma_I_ci_vec[2]) && sigma_I_ci_vec[1] > sigma_I_ci_vec[2]) {
          std::swap(sigma_I_ci_vec[1], sigma_I_ci_vec[2]);
        }
      }
    }
    
    // Confidence interval for sigma_G (sigma_a - top level)
    if (!ISNAN(w2U) && w2U != 0 && !ISNAN(c2) && !ISNAN(c3) &&
        !ISNAN(G1) && !ISNAN(H2) && !ISNAN(G12) && 
        !ISNAN(H1) && !ISNAN(G2) && !ISNAN(H12)) {
      double V_L_G_base = std::pow(G1, 2.0) * std::pow(S1_squared, 2.0) + std::pow(H2, 2.0) * std::pow(c2, 2.0) * std::pow(S2_squared, 2.0) + G12 * c2 * S1_squared * S2_squared;
      double V_U_G_base = std::pow(H1, 2.0) * std::pow(S1_squared, 2.0) + std::pow(G2, 2.0) * std::pow(c2, 2.0) * std::pow(S2_squared, 2.0) + H12 * c2 * S1_squared * S2_squared;
      double V_L_G = NA_REAL, V_U_G = NA_REAL;
      
      if (c3 >= 0) {
        if(!ISNAN(G3)&&!ISNAN(G32)&&!ISNAN(G13_star)&&!ISNAN(H3)&&!ISNAN(H32)&&!ISNAN(H13)){
          V_L_G = V_L_G_base + std::pow(G3, 2.0) * std::pow(c3, 2.0) * std::pow(S3_squared, 2.0) + G32 * c3 * c2 * S3_squared * S2_squared + G13_star * c3 * S1_squared * S3_squared;
          V_U_G = V_U_G_base + std::pow(H3, 2.0) * std::pow(c3, 2.0) * std::pow(S3_squared, 2.0) + H32 * c3 * c2 * S3_squared * S2_squared + H13 * c3 * S1_squared * S3_squared;
        }
      } else { // c3 < 0
        if(!ISNAN(H3)&&!ISNAN(G13)&&!ISNAN(H23_star)&&!ISNAN(G3)&&!ISNAN(H13)){
          V_L_G = V_L_G_base + std::pow(H3, 2.0) * std::pow(c3, 2.0) * std::pow(S3_squared, 2.0) + G13 * std::abs(c3) * S1_squared * S3_squared + H23_star * std::abs(c3) * c2 * S2_squared * S3_squared;
          V_U_G = V_U_G_base + std::pow(G3, 2.0) * std::pow(c3, 2.0) * std::pow(S3_squared, 2.0) + H13 * std::abs(c3) * S1_squared * S3_squared + H23_star * std::abs(c3) * c2 * S2_squared * S3_squared;
        }
      }
      
      if (!ISNAN(V_L_G) && !ISNAN(V_U_G) && V_L_G >= 0 && V_U_G >= 0) {
        double term_G_var_hat = (S1_squared - c2 * S2_squared + c3 * S3_squared) / w2U;
        double lower_var_G = term_G_var_hat - std::sqrt(V_L_G) / w2U;
        double upper_var_G = term_G_var_hat + std::sqrt(V_U_G) / w2U;
        
        // Truncate sigma_G at zero
        sigma_G_ci_vec[1] = std::sqrt(std::max(0.0, lower_var_G));
        sigma_G_ci_vec[2] = std::sqrt(std::max(0.0, upper_var_G));
        if (!ISNAN(sigma_G_ci_vec[1]) && !ISNAN(sigma_G_ci_vec[2]) && sigma_G_ci_vec[1] > sigma_G_ci_vec[2]) {
          std::swap(sigma_G_ci_vec[1], sigma_G_ci_vec[2]);
        }
      }
    }
  }
  
  // Apply mult and grand_mean for CV if necessary
  if (output_type == "cv_ci") {
    if (ISNAN(grand_mean) || grand_mean == 0) {
      Rcpp::warning("Grand mean is NA or zero; CV CIs cannot be calculated and will be NA.");
      for(int k=0; k<3; ++k) { sigma_G_ci_vec[k]=NA_REAL; sigma_I_ci_vec[k]=NA_REAL; sigma_A_ci_vec[k]=NA_REAL; }
      for (int i = 0; i < sigma_i_persubj_point.size(); ++i) { sigma_i_persubj_point[i] = NA_REAL; sigma_a_persubj_point[i] = NA_REAL; }
    } else {
      for(int k=0; k<3; ++k) {
        if(!ISNAN(sigma_G_ci_vec[k])) sigma_G_ci_vec[k] = sigma_G_ci_vec[k] / grand_mean * mult;
        if(!ISNAN(sigma_I_ci_vec[k])) sigma_I_ci_vec[k] = sigma_I_ci_vec[k] / grand_mean * mult;
        if(!ISNAN(sigma_A_ci_vec[k])) sigma_A_ci_vec[k] = sigma_A_ci_vec[k] / grand_mean * mult;
      }
      for (int i = 0; i < sigma_i_persubj_point.size(); ++i) {
        if(!ISNAN(sigma_i_persubj_point[i])) sigma_i_persubj_point[i] = sigma_i_persubj_point[i] / grand_mean * mult; else sigma_i_persubj_point[i] = NA_REAL;
        if(!ISNAN(sigma_a_persubj_point[i])) sigma_a_persubj_point[i] = sigma_a_persubj_point[i] / grand_mean * mult; else sigma_a_persubj_point[i] = NA_REAL;
      }
    }
  } else { // "sigma_ci"
    for(int k=0; k<3; ++k) {
      if(!ISNAN(sigma_G_ci_vec[k])) sigma_G_ci_vec[k] *= mult;
      if(!ISNAN(sigma_I_ci_vec[k])) sigma_I_ci_vec[k] *= mult;
      if(!ISNAN(sigma_A_ci_vec[k])) sigma_A_ci_vec[k] *= mult;
    }
    for (int i = 0; i < sigma_i_persubj_point.size(); ++i) {
      if(!ISNAN(sigma_i_persubj_point[i])) sigma_i_persubj_point[i] *= mult;
      if(!ISNAN(sigma_a_persubj_point[i])) sigma_a_persubj_point[i] *= mult;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("sigma_A") = sigma_A_ci_vec,
                            Rcpp::Named("sigma_I") = sigma_I_ci_vec,
                            Rcpp::Named("sigma_G") = sigma_G_ci_vec,
                            Rcpp::Named("sigma_i") = sigma_i_persubj_point,
                            Rcpp::Named("sigma_a") = sigma_a_persubj_point,
                            Rcpp::Named("HBHR") = HBHR,
                            Rcpp::Named("beta") = grand_mean);
}

//' @title Bootstrap Confidence Intervals for ANOVA Variance Components
//'
//' @name bv_anova_bootstrap_ci
//' 
//' @param data_orig A \code{list} or \code{data.frame} containing the original data.
//'   It must include the columns \code{SubjectID} (as integer), \code{SampleID}
//'   (as integer), \code{ReplicateID} (as integer), and \code{y} (as numeric).
//' @param B An \code{integer} specifying the number of bootstrap replicates to perform.
//'   A larger number of replicates (e.g., 1000 or more) is recommended for stable
//'   confidence intervals.
//' @param level A \code{double} representing the desired confidence level for the
//'   intervals (e.g., 0.95 for 95\% CIs). Defaults to 0.95.
//' @param output_type_for_point_est A \code{character} string specifying the metric
//'   to be estimated in each bootstrap iteration. Must be either \code{"sigma"}
//'   (for standard deviations) or \code{"cv"} (for coefficients of variation).
//'   The resulting CIs will correspond to this metric. Defaults to \code{"sigma"}.
//' @param mult A \code{double} used as a scaling factor for the estimates, applied
//'   in each bootstrap iteration. For example, set to 100 to express CVs as
//'   percentages. Defaults to 1.0.
//'   
//' @description
//' Calculates percentile bootstrap confidence intervals for variance components
//' estimated from a three-level nested ANOVA model.
//'
//' @details
//' This function implements a non-parametric case-resampling bootstrap
//' procedure. The key steps are:
//' \enumerate{
//'   \item Subjects (the highest level of hierarchy) are resampled with
//'         replacement from the original dataset.
//'   \item For each of the \code{B} bootstrap replicates, a new dataset is
//'         constructed by pooling all observations from the sampled subjects.
//'   \item The \code{\link{variance_components}} function is called on this
//'         new dataset to obtain point estimates of the variance components
//'         (\code{sigma_A}, \code{sigma_I}, \code{sigma_G}) and the HBHR.
//'   \item After all \code{B} replicates are completed, the confidence
//'         intervals are calculated from the percentiles of the resulting
//'         bootstrap distributions.
//' }
//' This case-resampling approach correctly accounts for the correlation
//' structure within subjects. The process can be interrupted by the user.
//' If any bootstrap iteration fails (e.g., due to a degenerate sample),
//' it is skipped, and a warning is suppressed to avoid cluttering the console.
//'
//' @return
//' A \code{list} with the following structure:
//' \describe{
//'   \item{\code{point_estimates}}{A \code{list} containing the point estimates of \code{sigma_A},
//'     \code{sigma_I}, \code{sigma_G}, and \code{HBHR} calculated from the original, un-resampled data.}
//'   \item{\code{conf_intervals}}{A \code{list} containing the lower and upper bounds of the
//'     percentile bootstrap confidence intervals for each component.}
//'   \item{\code{bootstrap_replicates}}{A \code{list} containing the raw numeric vectors of all
//'     \code{B} estimates for each component. This allows for manual inspection of the
//'     bootstrap distributions (e.g., plotting histograms).}
//'   \item{\code{n_bootstrap_replicates}}{The number of bootstrap replicates (\code{B}) performed.}
//'   \item{\code{confidence_level}}{The confidence level used for the intervals.}
//' }
//'
//' @examples
//' # Create a sample data list for demonstration
//' sample_data_list <- list(
//'   SubjectID = as.integer(rep(1:3, each = 4)),
//'   SampleID = as.integer(rep(1:6, each = 2)),
//'   ReplicateID = as.integer(rep(1:2, 6)),
//'   y = c(
//'     rnorm(2, 100, 10), rnorm(2, 105, 10), # Subject 1
//'     rnorm(2, 120, 12), rnorm(2, 118, 12), # Subject 2
//'     rnorm(2, 90, 8), rnorm(2, 92, 8)      # Subject 3
//'   )
//' )
//'
//' # Run bootstrap with a small number of replicates for the example
//' # In practice, B should be >= 1000
//' # Note: This example requires the 'variance_components' function to be available.
//' # The call is commented out as it can be slow and depends on other functions.
//'
//' # bootstrap_results <- bv_anova_bootstrap_ci(
//' #   data_orig = sample_data_list,
//' #   B = 50, # Use a small B for speed in an example
//' #   level = 0.95,
//' #   output_type_for_point_est = "sigma"
//' # )
//'
//' # print(bootstrap_results$point_estimates)
//' # print(bootstrap_results$conf_intervals)
//'

// [[Rcpp::export]]
List bv_anova_bootstrap_ci(
    List data_orig,
    int B,
    double level = 0.95,
    std::string output_type_for_point_est = "sigma",
    double mult = 1.0
) {
  
  // --- 1. Extract original data and identify unique subjects ---
  IntegerVector subjects_orig_all = data_orig["SubjectID"];
  IntegerVector samples_orig_all = data_orig["SampleID"];
  IntegerVector replicates_orig_all = data_orig["ReplicateID"];
  NumericVector values_orig_all = data_orig["y"];
  
  IntegerVector unique_subject_ids_vec = Rcpp::unique(subjects_orig_all).sort();
  int n_unique_subjects = unique_subject_ids_vec.size();
  
  // --- 2. Storage for bootstrap estimates ---
  NumericMatrix boot_estimates_sigma_A(B, 1);
  NumericMatrix boot_estimates_sigma_I(B, 1);
  NumericMatrix boot_estimates_sigma_G(B, 1);
  NumericMatrix boot_estimates_HBHR(B, 1);
  NumericMatrix boot_estimates_beta(B, 1);
  
  
  // --- 3. Bootstrap Loop ---
  for (int b = 0; b < B; ++b) {
    Rcpp::checkUserInterrupt(); // Allow user to interrupt
    
    // --- 3a. Resample Subject IDs (indices of unique_subject_ids_vec) ---
    IntegerVector sampled_indices(n_unique_subjects);
    for(int i=0; i < n_unique_subjects; ++i) {
      sampled_indices[i] = R::runif(0, n_unique_subjects); // Sample index with replacement
    }
    
    // --- 3b. Construct Bootstrap Dataset ---
    IntegerVector boot_subj_id_col;
    IntegerVector boot_samp_id_col;
    IntegerVector boot_repl_id_col;
    NumericVector boot_y_col;
    
    for (int i = 0; i < n_unique_subjects; ++i) {
      int original_subj_id_to_include = unique_subject_ids_vec[sampled_indices[i]];
      int new_bootstrap_subj_id = i + 1; // Assign new contiguous IDs for this bootstrap sample
      
      for (int k = 0; k < subjects_orig_all.size(); ++k) {
        if (subjects_orig_all[k] == original_subj_id_to_include) {
          boot_subj_id_col.push_back(new_bootstrap_subj_id);
          boot_samp_id_col.push_back(samples_orig_all[k]);
          boot_repl_id_col.push_back(replicates_orig_all[k]);
          boot_y_col.push_back(values_orig_all[k]);
        }
      }
    }
    
    if (boot_subj_id_col.size() == 0) { // Should not happen if n_unique_subjects > 0
      boot_estimates_sigma_A(b, 0) = NA_REAL;
      boot_estimates_sigma_I(b, 0) = NA_REAL;
      boot_estimates_sigma_G(b, 0) = NA_REAL;
      boot_estimates_HBHR(b, 0) = NA_REAL;
      boot_estimates_beta(b, 0) = NA_REAL;
      continue;
    }
    
    List data_boot = List::create(
      Named("SubjectID") = boot_subj_id_col,
      Named("SampleID") = boot_samp_id_col,
      Named("ReplicateID") = boot_repl_id_col,
      Named("y") = boot_y_col
    );
    
    // --- 3c. Estimate Parameters on Bootstrap Dataset ---
    try {
      // We need point estimates here, so we use "sigma" or "cv" output type
      // The 'level' argument to variance_components is ignored for "sigma"/"cv"
      List vc_results_boot = variance_components(data_boot, output_type_for_point_est, mult, 0.95);
      
      // The variance_components function returns scalar for "sigma" or "cv" output
      boot_estimates_sigma_A(b, 0) = as<double>(vc_results_boot["sigma_A"]);
      boot_estimates_sigma_I(b, 0) = as<double>(vc_results_boot["sigma_I"]);
      boot_estimates_sigma_G(b, 0) = as<double>(vc_results_boot["sigma_G"]);
      boot_estimates_HBHR(b, 0) = as<double>(vc_results_boot["HBHR"]);
      boot_estimates_beta(b, 0) = as<double>(vc_results_boot["beta"]);
      
      
    } catch (std::exception &ex) {
      // If bv_anova or variance_components fails (e.g., due to all NA or singular data in a bootstrap sample)
      boot_estimates_sigma_A(b, 0) = NA_REAL;
      boot_estimates_sigma_I(b, 0) = NA_REAL;
      boot_estimates_sigma_G(b, 0) = NA_REAL;
      boot_estimates_HBHR(b, 0) = NA_REAL;
      boot_estimates_beta(b, 0) = NA_REAL;
      // Rf_warning("Bootstrap iteration %d failed: %s", b + 1, ex.what());
    } catch (...) {
      boot_estimates_sigma_A(b, 0) = NA_REAL;
      boot_estimates_sigma_I(b, 0) = NA_REAL;
      boot_estimates_sigma_G(b, 0) = NA_REAL;
      boot_estimates_HBHR(b, 0) = NA_REAL;
      boot_estimates_beta(b, 0) = NA_REAL;
      // Rf_warning("Bootstrap iteration %d failed due to unknown C++ error.", b + 1);
    }
  }
  
  // --- 4. Calculate Percentile Confidence Intervals ---
  double alpha = 1.0 - level;
  double lower_prob = alpha / 2.0;
  double upper_prob = 1.0 - alpha / 2.0;
  
  Function Rquantile("quantile"); // Get R's quantile function
  
  NumericVector sigma_A_ci_res(2, NA_REAL);
  NumericVector sigma_I_ci_res(2, NA_REAL);
  NumericVector sigma_G_ci_res(2, NA_REAL);
  NumericVector HBHR_ci_res(2, NA_REAL);
  NumericVector beta_ci_res(2, NA_REAL);
  
  NumericVector probs = NumericVector::create(lower_prob, upper_prob);
  
  // Helper lambda to compute CIs
  auto get_ci = [&](NumericVector estimates) {
    NumericVector valid_estimates;
    for(int i=0; i<estimates.size(); ++i) {
      if(!R_IsNA(estimates[i]) && R_finite(estimates[i])) { // Check for NA and Inf/-Inf
        valid_estimates.push_back(estimates[i]);
      }
    }
    if (valid_estimates.size() < 2) { // Need at least 2 points for percentile
      return NumericVector::create(NA_REAL, NA_REAL);
    }
    return as<NumericVector>(Rquantile(valid_estimates, probs, Named("na.rm", true), Named("type", 7))); // type 7 is R default
  };
  
  sigma_A_ci_res = get_ci(boot_estimates_sigma_A(_,0)); // Extract column
  sigma_I_ci_res = get_ci(boot_estimates_sigma_I(_,0));
  sigma_G_ci_res = get_ci(boot_estimates_sigma_G(_,0));
  HBHR_ci_res = get_ci(boot_estimates_HBHR(_,0));
  beta_ci_res = get_ci(boot_estimates_beta(_,0));
  
  // Also return point estimates from original data for convenience
  List vc_results_orig = variance_components(data_orig, output_type_for_point_est, mult, level);
  
  return List::create(
    Named("point_estimates") = List::create(
      Named("sigma_A") = as<double>(vc_results_orig["sigma_A"]),
      Named("sigma_I") = as<double>(vc_results_orig["sigma_I"]),
      Named("sigma_G") = as<double>(vc_results_orig["sigma_G"]),
      Named("HBHR") = as<double>(vc_results_orig["HBHR"]),
      Named("beta") = as<double>(vc_results_orig["beta"])
    ),
    Named("conf_intervals") = List::create(
      Named("sigma_A_CI") = sigma_A_ci_res,
      Named("sigma_I_CI") = sigma_I_ci_res,
      Named("sigma_G_CI") = sigma_G_ci_res,
      Named("HBHR_CI")    = HBHR_ci_res,
      Named("beta_CI")    = beta_ci_res
    )
  );
}