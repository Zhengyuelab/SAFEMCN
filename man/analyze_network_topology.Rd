\name{analyze_network_topology}
\alias{analyze_network_topology}
\title{Analyze Network Topology Parameters with Rarefaction}
\description{
Calculate network topology parameters from OTU tables with customizable correlation thresholds, parallel processing options, and visualization capabilities including topology parameter plots and exponential formula fitting.
}
\usage{
analyze_network_topology(work_dir, otu_file, r_threshold = 0.6, 
                        p_threshold = 0.05, cor_method = "spearman",
                        start_size = 5, end_size = NULL, step_size = 1,
                        replicates = 50, use_parallel = FALSE,
                        output_prefix = "network_parameters",
                        plot_topology = FALSE, plot_start_size = 11,
                        fit_formula = FALSE, fit_start_size = 11)
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
  \item{fit_formula}{Whether to perform exponential formula fitting (default: FALSE)}
  \item{fit_start_size}{Starting sample size for formula fitting (default: 11)}
}
\value{
A list containing raw results and summary statistics
}
\examples{
results <- analyze_network_topology(
  work_dir = "/path/to/directory",
  otu_file = "otu_table.csv",
  r_threshold = 0.6,
  p_threshold = 0.05,
  cor_method = "spearman",
  start_size = 5,
  end_size = 58,
  step_size = 1,
  replicates = 50,
  use_parallel = TRUE,
  output_prefix = "my_network",
  plot_topology = TRUE,
  plot_start_size = 11,
  fit_formula = TRUE,
  fit_start_size = 11
)
}