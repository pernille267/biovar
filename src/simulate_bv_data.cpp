#include <Rcpp.h>
using namespace Rcpp;

//' Simulate Biological Variation Data
//'
//' @title Simulate Biological Variation Data
//' @name simulate_bv_data
//'
//' @param n \code{Integer}. The number of subjects.
//' @param S \code{Integer}. The number of samples taken from each subject. Defaults to 10.
//' @param R \code{Integer}. The number of replicates measured for each sample. Defaults to 2.
//' @param cvi \code{Double}. The coefficient of variation (in percent) for variability within each subject (default is 10).
//' @param cva \code{Double}. The coefficient of variation (in percent) for the variability of analytical error (default is 2). This is the source of variation between replicated measurements.
//' @param cvg \code{Double}. The coefficient of variation (in percent) for variability between subjects (default is 50).
//' @param mu \code{Double}. The overall true mean for the population (default is 100).
//'
//' @description This function simulates biological variation data based on a specified number of subjects, samples, and replicates.
//' @details This function generates simulated biological variation data. It first generates subject-specific means based on a normal distribution with mean `mu` and coefficient of variation `cvg`. For each subject, it then generates `S` samples based on their specific mean and the coefficient of variation `cvi`. Finally, for each sample, it generates `R` replicate measurements based on the sample mean and the coefficient of variation `cva`.
//'
//' @return A `list` with columns:
//' \itemize{
//'   \item \code{SubjectID} - The ID of the subject.
//'   \item \code{SampleID} - The ID of the sample.
//'   \item \code{ReplicateID} - The ID of the replicate.
//'   \item \code{y} - The measured value for the replicate.
//' }
//'
//' @examples
//' \dontrun{
//' // Simulate data for 100 subjects, with 10 samples per subject and 3 replicates per sample
//' output <- simulate_bv_data(15, 10, 2, 10, 2, 50, 100)
//' output <- as.data.table(output)
//' print(output)
//' }
//'

// [[Rcpp::export]]
List simulate_bv_data(int n, int S = 10, int R = 2, double cvi = 10, double cva = 2, double cvg = 50, double mu) {
  // Set up vectors to store the results
  std::vector<int> SubjectID(n * S * R);
  std::vector<int> SampleID(n * S * R);
  std::vector<int> ReplicateID(n * S * R);
  std::vector<double> y(n * S * R);
  
  // Generate subjects
  std::vector<double> subjects(n);
  for(int i = 0; i < n; i++) {
    subjects[i] = std::abs(R::rnorm(mu, (cvg / 100.0) * mu));
  }
  
  // Generate samples and replicates
  int idx = 0;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < S; j++) {
      double sample = std::abs(R::rnorm(subjects[i], (cvi / 100.0) * subjects[i]));
      for(int k = 0; k < R; k++) {
        double replicate = std::abs(R::rnorm(sample, (cva / 100.0) * sample));
        SubjectID[idx] = i + 1;
        SampleID[idx] = j + 1;
        ReplicateID[idx] = k + 1;
        y[idx] = replicate;
        idx++;
      }
    }
  }
  
  // Create a list to return
  return List::create(Named("SubjectID") = SubjectID,
                      Named("SampleID") = SampleID,
                      Named("ReplicateID") = ReplicateID,
                      Named("y") = y);
}

