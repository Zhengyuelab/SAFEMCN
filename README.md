# SAFEMCN <img src="https://raw.githubusercontent.com/Zhengyuelab/SAFEMCN/main/logo2.png" width="80" align="right" />

SAFEMCN:An R package implementing SAFEMCN, a sample-size-aware framework that integrates rarefaction-style resampling with lag-1 autocorrelation (AR1) diagnostics to evaluate the reliability and reproducibility of microbial co-occurrence networks across sample-size gradients.


# **Overview**  
SAFEMCN is a comprehensive R language software package designed to address the critical issue of sample size sufficiency in microbial ecology. It offers a reliable and statistically-based process for generating network topological parameters under different sample sizes through random sampling, and uses exponential fitting to describe the dynamic changes in network properties under the gradient of sample size. It can also use the fitting formula to calculate the network topological parameters for larger sample sizes. On this basis, the minimum effective sample size (Nmin) is determined through stability analysis based on AR1, ensuring the construction of a repeatable and effective microbial co-occurrence network.

<p align="center">
  <img src="https://raw.githubusercontent.com/Zhengyuelab/SAFEMCN/main/logo.png" 
       width="65%" 
       alt="Workflow of the SAFEMCN R package" />
</p>
<p align="center">
  <b>Figure 1. Workflow of the SAFEMCN R package</b>
</p>

# **Main Features**  
* From ASV quality control to network inference.  
* Customizable correlation methods (Spearman/Pearson) and significance thresholds (R and p value).  
* Iterative bootstrapping to track topological changes across sample-size gradients.  
* The exponential fitting method is employed to fit the variation of network topology parameters with the sample size.  
* This fitting formula is used to predict the network topology parameters under large sample sizes.
* The lag-1 autocorrelation (AR1) is calculated to determine the minimum effective sample size (Nmin) for network structure stability.

# **Install R/RStudio**  
If you do not already have R/RStudio installed, follow these steps:   
1. Install [R](https://www.r-project.org/)
2. Install [RStudio](https://posit.co/download/rstudio-desktop/)

Open RStudio -> Tools -> Global Options -> Packages, select the appropriate mirror in Primary CRAN repository.

# **Install SAFEMCN** 
Install SAFEMCN package from CRAN.
```r
install.packages("SAFEMCN")
```
Install the latest development version from Github.
```r
#If devtools package is not installed, first install it.
install.packages("devtools")
devtools::install_github("Zhengyuelab/SAFEMCN")
```

# **Usage**  
### **1. Code Execution and Parameter Settings**
```r
library(SAFEMCN)

results <- analyze_network_topology(
  work_dir = "./",
  otu_file = "otu.csv",
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
  ar1_windows = c(10, 15)
)
```
* work_dir: Specifies the working directory for reading input files and saving output results.  
* otu_file: Input OTU/ASV table (rows: features, cols: samples).  
* r_threshold: The correlation coefficient threshold (|R|).   
* p_threshold: The significance level(filters non-significant edges). 
* cor_method: Correlation algorithm ("spearman" or "pearson").  
* start_size & end_size: Defines the range of the sample-size gradient (e.g., from 5 to 58 samples).  
* step_size: Increment of gradient (1 = single sample steps).  
* replicates: The number of repeated sampling for each sample size (e.g., 50 repeats).
* use_parallel: Enable multi-core parallel processing.
* output_prefix: Prefix for exported result files. 
* plot_topology: Logical; if TRUE, generates plots showing the trend of topological parameters across the sample-size gradient.  
* plot_start_size: Start sample size for plotting.  
* fit_formula: Logical; if TRUE, Apply exponential fitting to parameters.  
* fit_start_size:  Start sample size for fitting.  
* predict_end_size: Target sample size for extrapolation.  
* run_ar1: Logical; if TRUE, performs lag-1 Autoregressive (AR1) analysis to calculate the AR1 coefficient at each sample size, serving as an indicator of network stability.  
* ar1_input_file: The input file for AR1 analysis, containing network topological parameters (e.g., node number, edge number, or density) calculated across the sample-size gradient.  
* ar1_windows: Defines the sliding window sizes (e.g., 10 and 15) for calculating rolling statistics like standard deviation or autocorrelation in the AR1 stability analysis.

### **2. Case Studies**
#### **2.1 Aloha Dataset Example**
Based on the long-term time-series dataset from Station ALOHA, a gradient of microbial co-occurrence networks was constructed by varying the sample size from 11 to 58 with a step size of 1. As illustrated in Figure 2, key network topological features—including network nodes, edges, density, diameter, average path length, and average degree—exhibited distinct dynamic patterns with expanding sample sizes, based on which the minimum sample size (N<sub>min</sub>) required for reliable network construction was successfully determined via lag-1 autoregressive (AR1) coefficients.

<p align="center">
  <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/ALOHA.png" 
       width="65%" 
       alt="Variations in network topological parameters and the lag-1 autocorrelation (AR1) coefficient" />
</p>
<p align="center">
  <b>Figure 2. Variations in network topological parameters and the lag-1 autocorrelation (AR1) coefficient</b>
</p>

#### 📌 **Result Interpretation:**
With increasing sample sizes, node numbers expanded toward saturation, whereas edge counts dramatically decreased before leveling off. Across these trajectories, the network topology parameters stabilized in a clear sequential order: node first, followed by edge, density, average degree, and lastly, diameter and average path length, reflecting a gradient of sensitivity to sample size.
The Nmin threshold was further identified by the AR1 module of the SAFEMCN package. The rapid decrease and then rapid increase in the AR1 value indicated the emergence of the topological plateau, providing a basis to determine the critical point corresponding to N<sub>min</sub>. To provide a robust estimate, the final Nmin was calculated as the average of the N<sub>min</sub> derived from all six topological parameters, resulting in a value of 31 ± 5. This threshold marked the earlier sampling depth at which inferred topology became reproducible across resampling replicates and could be regarded as meeting a minimum level of statistical validity.

#### **2.2 Exponential Fitting and Prediction**


<p align="center" style="margin-bottom: 2px;">
  <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/Predicting_Nmin.png" 
       width="55%" 
       alt="Sample-size pre-estimation for microbial co-occurrence networks" />

<p align="center" style="margin-top: 0; margin-bottom: 5px;">
  <b>Figure 3. Sample-size pre-estimation for microbial co-occurrence networks</b>

