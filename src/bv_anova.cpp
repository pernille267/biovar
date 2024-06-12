#include <Rcpp.h>
using namespace Rcpp;

// Taking mean and remove NA values
double mean_no_na(NumericVector x) {
  double sum = 0;
  int count = 0;
  int n = x.size();
  if(n == 1){
    return x[0];
  }
  
  for(int i = 0; i < x.size(); i++) {
    if(!NumericVector::is_na(x[i])) {
      sum += x[i];
      count++;
    }
  }
  return count == 0 ? NA_REAL : sum / count;
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
//' 
//'
//' @description This function uses biological variation data to estimate both one-way (one model for each subject) and two-way (general model overall subjects) nested ANOVA models. 
//' @details The components of the estimated model can be used to estimate intra biological variation, inter biological variation and analytical variation, with corresponding confidence intervals. The model estimates relies on that all model effects are mutually independent and normally distributed.
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
List bv_anova(List data) {
  IntegerVector subjects = data["SubjectID"];
  IntegerVector samples = data["SampleID"];
  IntegerVector replicates = data["ReplicateID"];
  NumericVector values = data["y"];
  int N = subjects.size();
  int n = max(subjects);
  double grand_mean = mean_no_na(values);
  double weighted_grand_mean = 0;
  bool grand_mean_is_na = ISNAN(grand_mean);
  if(grand_mean_is_na){
    Rcpp::stop("Grand mean is missing...");
  }
  
  double SST = 0;
  double SS1 = 0;
  double SS2 = 0;
  double SS1U = 0;
  double SS2U = 0;
  double w1U = 0;
  double w1U_num = 0;
  double w1U_denom = 0;
  double w2U = 0;
  double w2U_denom = 0;
  double w3U = 0;
  double w3U_denom = 0;
  int nT = -1;
  int n1 = n - 1;
  int n2 = -n;
  
  
  NumericVector weighted_subject_means(n);
  NumericVector K_i(n);
  NumericVector SSi(n);
  NumericVector SSt(n);
  NumericVector SSa(n);
  IntegerVector ni(n);
  IntegerVector na(n);
  NumericVector Si_squared(n);
  NumericVector Sa_squared(n);
  
  for (int k = 0; k < N; ++k) {
    SST += pow(values[k] - grand_mean, 2);
  }
  
  for (int i = 1; i <= n; ++i) {
    double R_mean_i = 0;
    double R_mean_denom = 0;
    double weighted_subject_mean_i = 0;
    
    // Extract relevant information from subject i
    IntegerVector subject_i_samples = subset_samples(samples, subjects, i);
    NumericVector subject_i_values = subset_samples_values(values, subjects, i);
    int S_i = max(subject_i_samples);
    n2 += S_i;
    NumericVector weighted_sample_means(S_i);
    
    // Initialize degrees of freedom for variability within subject i
    ni[i - 1] = S_i - 1;
    na[i - 1] = -S_i;
    
    // Calculate mean of all values within subject i
    double subject_mean_i = mean_no_na(subject_i_values);
    bool subject_mean_is_na = ISNAN(subject_mean_i);
    if(subject_mean_is_na){
      Rcpp::stop("Subject mean is missing...");
    }
    
    // Calculate SST within subject i
    for (int k = 0; k < subject_i_samples.size(); ++k) {
      SSt[i - 1] += pow(subject_i_values[k] - subject_mean_i, 2);
    }
    
    for (int j = 1; j <= S_i; ++j) {
      
      // Extract relevant information from sample j within subject i
      IntegerVector subject_i_sample_j_replicates = subset_replicates(replicates, subjects, samples, i, j);
      NumericVector subject_i_sample_j_values = subset_replicates_values(values, subjects, samples, i, j);
      double R_ij = max(subject_i_sample_j_replicates);
      nT += R_ij;
      na[i - 1] += R_ij;
      R_mean_denom += 1.0 / R_ij;
      SS1 += R_ij * pow(subject_mean_i - grand_mean, 2);
      
      // Calculate mean of all values within sample j within subject i
      double subject_i_mean_sample_j = mean_no_na(subject_i_sample_j_values);
      bool subject_i_mean_sample_j_is_na = ISNAN(subject_i_mean_sample_j);
      if(subject_i_mean_sample_j_is_na){
        Rcpp::stop("Subject i sample j mean is missing...");
      }
      
      weighted_sample_means[j - 1] = subject_i_mean_sample_j;
      weighted_subject_mean_i += subject_i_mean_sample_j / S_i;
      
      SS2 += R_ij * pow(subject_i_mean_sample_j - subject_mean_i, 2);
      SSi[i - 1] += R_ij * pow(subject_i_mean_sample_j - subject_mean_i, 2);
    }
    
    SSa[i - 1] = SSt[i - 1] - SSi[i - 1];
    Si_squared[i - 1] = SSi[i - 1] / ni[i - 1];
    Sa_squared[i - 1] = SSa[i - 1] / na[i - 1];
    
    weighted_grand_mean += weighted_subject_mean_i / n;
    weighted_subject_means[i - 1] = weighted_subject_mean_i;
    
    R_mean_i = S_i / R_mean_denom;
    w2U_denom += 1.0 / S_i / R_mean_i;
    w3U_denom += (S_i - 1.0) / R_mean_i;
    w1U_num += 1.0 / S_i;
    w1U_denom += 1.0 / S_i / R_mean_i;
    
    for (int j = 0; j < S_i; ++j){
      SS2U += pow(weighted_sample_means[j] - weighted_subject_mean_i, 2);
    }
  }
  w1U = w1U_num / w1U_denom;
  w2U = n / w2U_denom;
  w3U = n2 / w3U_denom;
  
    
  for(int i = 0; i < n; ++i){
    SS1U += pow(weighted_subject_means[i] - weighted_grand_mean, 2);
  }
  
  SS1U = SS1U * w2U;
  SS2U = SS2U * w3U;
  
  double SS3 = SST - SS1 - SS2;
  int n3 = nT - n1 - n2;
  
  double S1_squared = SS1 / n1;
  double S2_squared = SS2 / n2;
  double S3_squared = SS3 / n3;
  double S1U_squared = SS1U / n1;
  double S2U_squared = SS2U / n2;
  
  return List::create(Named("S1_squared") = S1_squared,
                      Named("S2_squared") = S2_squared,
                      Named("S3_squared") = S3_squared,
                      Named("S1U_squared") = S1U_squared,
                      Named("S2U_squared") = S2U_squared,
                      Named("n1") = n1,
                      Named("n2") = n2,
                      Named("n3") = n3,
                      Named("w1U") = w1U,
                      Named("w2U") = w2U,
                      Named("w3U") = w3U,
                      Named("Si_squared") = Si_squared,
                      Named("Sa_squared") = Sa_squared,
                      Named("ni") = ni,
                      Named("na") = na);
}
