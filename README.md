# SAFEMCN
SAFEMCN:An R package for evaluating the reliability and reproducibility of microbial co-occurrence networks across sample-size gradients.

# **Overview**  
SAFEMCN is a comprehensive R package designed to address the critical challenge of sample size adequacy in microbial ecology. It provides a robust, statistically grounded workflow to generate topological distributions through resampling and characterizes the dynamic changes of network properties across sample-size gradients using exponential fitting.  
On this basis, the minimum effective sample size (Nmin) is determined by subsequent AR1-based stability analysis implemented in Python, ensuring the construction of reproducible and valid microbial co-occurrence networks.

# **Main Features**  
From ASV quality control to network inference.  
Customizable correlation methods (Spearman/Pearson) and significance thresholds (R and p value).  
Iterative bootstrapping to track topological changes across sample-size gradients.  
The exponential fitting method is employed to fit the variation of network topology parameters with the sample size.  

# **Install R/RStudio**  
If you do not already have R/RStudio installed, follow these steps:   
Install R  
Install RStudio  
Open RStudio -> Tools -> Global Options -> Packages, select the appropriate mirror in Primary CRAN repository.

# **Install SAFEMCN**  
Install the latest development version from Github.  
If devtools package is not installed, first install it  
```r
install.packages("devtools")
```
```r
devtools::install_github("Zhengyuelab/SAFEMCN")
```

# **Usage**  
```r
library(SAFEMCN)

results <- analyze_network_topology(  
  work_dir = "./",  
  otu_file = "otu.csv",  
  r_threshold = 0.6,  
  p_threshold = 0.05,  
  cor_method = "spearman",  
  start_size = 5,  
  end_size = 58,  
  step_size = 1,  
  replicates = 50,  
  use_parallel = TRUE,  
  output_prefix = "network_parameters",  
  plot_topology = TRUE,  
  plot_start_size = 11,  
  fit_formula = TRUE,
  fit_start_size = 11
)
```
work_dir: Specifies the working directory for reading input files and saving output results.  
otu_file: The input ASV/OTU abundance table (typically a .csv file where rows are features and columns are samples).  
r_threshold: The correlation coefficient threshold (|R|).   
p_threshold: The significance level. Used to filter out non-significant correlations to ensure edge reliability.  
cor_method: The algorithm for calculating correlations. Options include "spearman" (robust for non-normal distributions) or "pearson".  
start_size & end_size: Defines the range of the sample-size gradient (e.g., from 5 to 58 samples).  
step_size: The increment for the gradient. A step of 1 means the analysis is performed for every single sample increment.  
replicates: The number of bootstrapping iterations per sample size (e.g., 50 repeats) to generate a distribution of topological metrics and reduce stochastic error.  
use_parallel: Enables multi-core parallel processing to significantly accelerate the heavy computation required for thousands of network reconstructions.  
output_prefix: The prefix for exported files, including CSV data for subsequent Python-based AR1 analysis and visualization plots.  
plot_topology: Logical; if TRUE, generates plots showing the trend of topological parameters across the sample-size gradient.  
plot_start_size: The starting sample size for visualization, often used to exclude highly volatile results from very small sample sizes.  
fit_formula: Logical; if TRUE, applies exponential fitting to the data to characterize the mathematical convergence of network properties.  
fit_start_size: Specifies the sample size at which the fitting process begins to ensure a biologically meaningful and stable model.  
