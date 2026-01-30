# SAFEMCN
SAFEMCN:An R package for evaluating the reliability and reproducibility of microbial co-occurrence networks across sample-size gradients.

**Overview**  
SAFEMCN is a comprehensive R package designed to address the critical challenge of sample size adequacy in microbial ecology. It provides a robust, statistically grounded workflow to generate topological distributions through resampling and characterizes the dynamic changes of network properties across sample-size gradients using exponential fitting.
On this basis, the minimum effective sample size (Nmin) is determined by subsequent AR1-based stability analysis implemented in Python, ensuring the construction of reproducible and valid microbial co-occurrence networks.

**Main Features**  
From ASV quality control to network inference.
Customizable correlation methods (Spearman/Pearson) and significance thresholds (R and p value).
Iterative bootstrapping to track topological changes across sample-size gradients.
The exponential fitting method is employed to fit the variation of network topology parameters with the sample size.

**Install R/RStudio**  
If you do not already have R/RStudio installed, follow these steps: 
Install R
Install RStudio
Open RStudio -> Tools -> Global Options -> Packages, select the appropriate mirror in Primary CRAN repository.

**Install SAFEMCN**  
Install the latest development version from Github.
# If devtools package is not installed, first install it
install.packages("devtools")
devtools::install_github("Zhengyuelab/SAFEMCN")

**Usage**  
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
