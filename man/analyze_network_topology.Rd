\name{analyze_network_topology}
\alias{analyze_network_topology}
\title{Analyze Network Topology Parameters with Rarefaction}
\description{
Calculate network topology parameters from OTU tables with customizable correlation thresholds, parallel processing options, and visualization capabilities including topology parameter plots, exponential formula fitting (exp3P via nls), prediction of future sample sizes, combined/original CSV generation, and AR1 analysis.
}
\usage{
analyze_network_topology(work_dir, otu_file, r_threshold = 0.6,
                        p_threshold = 0.05, cor_method = "spearman",
                        start_size = 5, end_size = NULL, step_size = 1,
                        replicates = 50, use_parallel = FALSE,
                        output_prefix = "network_parameters",
                        plot_topology = FALSE, plot_start_size = 11,
                        fit_formula = FALSE, fit_start_size = 11,
                        predict_end_size = NULL,
                        run_ar1 = FALSE, ar1_input_file = NULL,
                        ar1_windows = c(5, 10, 15))
}
\arguments{
  \item{work_dir}{Working directory path}
  \item{otu_file}{OTU table file name}
  \item{r_threshold}{Correlation coefficient threshold (default: 0.6)}
  \item{p_threshold}{P-value threshold (default: 0.05)}
  \item{cor_method}{Correlation method: "spearman" or "pearson" (default: "spearman")}
  \item{start_size}{Starting sample size (default: 5)}
  \item{end_size}{Ending sample size (default: NULL, uses total samples)}
  \item{step_size}{Step size for sample sizes (default: 1)}
  \item{replicates}{Number of replicates (default: 50)}
  \item{use_parallel}{Use parallel processing (default: FALSE)}
  \item{output_prefix}{Output file prefix (default: "network_parameters")}
  \item{plot_topology}{Whether to plot topology parameters vs sample size (default: FALSE)}
  \item{plot_start_size}{Starting sample size for topology plots (default: 11)}
  \item{fit_formula}{Whether to perform exponential formula fitting with nls (default: FALSE)}
  \item{fit_start_size}{Starting sample size for formula fitting (default: 11)}
  \item{predict_end_size}{End sample size for prediction using fit formula (default: NULL, no prediction). When set, predicts topology parameters for sample sizes from end_size+1 to predict_end_size.}
  \item{run_ar1}{Whether to run AR1 (autocorrelation lag-1) analysis (default: FALSE)}
  \item{ar1_input_file}{Input file name for AR1 analysis (default: NULL, auto-selects node_num_combined.csv or node_num_original.csv). Can be any generated combined or original CSV file, e.g. "edge_num_combined.csv", "node_num_original.csv".}
  \item{ar1_windows}{Numeric vector of sliding window sizes for AR1 calculation (default: c(5, 10, 15))}
}
\value{
A list containing:
\item{raw_results}{Data frame of raw topology parameter results for all replicates}
\item{summary}{Data frame of mean and sd for each topology parameter at each sample size}
\item{fit_parameters}{(if fit_formula=TRUE) Data frame of exp3P fit parameters (a, b, c, R2) for each topology parameter}
\item{predictions}{(if predict_end_size is set) Data frame of predicted topology parameter values for future sample sizes}
}
\examples{
results <- analyze_network_topology(
  work_dir = "/path/to/directory",
  otu_file = "otu_table.csv",
  r_threshold = 0.6,
  p_threshold = 0.05,
  cor_method = "spearman",
  start_size = 5,
  end_size = 28,
  step_size = 1,
  replicates = 50,
  use_parallel = TRUE,
  output_prefix = "network_parameters",
  plot_topology = TRUE,
  plot_start_size = 11,
  fit_formula = TRUE,
  fit_start_size = 11,
  predict_end_size = 58,
  run_ar1 = TRUE,
  ar1_input_file = "node_num_combined.csv",
  ar1_windows = c(5, 10, 15)
)
}