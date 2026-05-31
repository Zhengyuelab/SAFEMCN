# *SAFEMCN* <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/logo.png" width="80" align="right" />

***SAFEMCN*** : *SAFEMCN* (A Sample-Size-Aware Framework for Evaluating Microbial Co-occurrence Networks) is dependent on R ≥ 3.5, the other dependencies in the R environment are igraph, Hmisc, parallel, ggplot2, grDevices, graphics, rlangstats and utils.


# **Overview**  
The *SAFEMCN* R package is a standardized framework designed to determine the minimum sample size (N<sub>min</sub>) required for robust microbial co-occurrence network inference, which systematically integrates three core functional modules: (1) Network construction through automated, parallelized subsampling across designated sample-size gradients; (2) Determination of the N<sub>min</sub> via lag-1 autoregressive (AR1) coefficient assessments; and (3) N<sub>min</sub> estimation via exponential model extrapolation, which provides a predictive, three-parameter exponential fitting strategy to prospectively forecast dataset requirements under constrained sample sizes.

<p align="center">
  <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/Workflow-SAFEMCN.png" 
       width="65%" 
       alt="Variations in network topological parameters and the lag-1 autocorrelation (AR1) coefficient" />
</p>


# **Main Features**  
* Traces topological structure changes across a sample-size gradient using a rarefaction-curve approach.
* The lag-1 autocorrelation (AR1) coefficient is calculated to determine the minimum effective sample size (N<sub>min</sub>) for network structure stability.
* Extrapolates topological trajectories beyond the sample-size range via a three-parameter exponential model and uses lag-1 autocorrelation (AR1) coefficients on fitted sequences to predict N<sub>min</sub>.
* Supports parallel computation for large-scale datasets by distributing iterations across multiple CPU cores.


# **Install R/RStudio**  
If you do not already have R/RStudio installed, follow these steps:   
1. Install [R](https://www.r-project.org/)
2. Install [RStudio](https://posit.co/download/rstudio-desktop/)

Open RStudio -> Tools -> Global Options -> Packages, select the appropriate mirror in Primary CRAN repository.

# **Install *SAFEMCN*** 
Install *SAFEMCN*  package from CRAN.
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
### **1. Parameter Settings**
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
* fit_formula: Logical; if TRUE, apply exponential fitting to parameters.  
* fit_start_size:  Start sample size for fitting.  
* predict_end_size: Target sample size for extrapolation.  
* run_ar1: Logical; if TRUE, performs lag-1 Autoregressive (AR1) analysis to calculate the AR1 coefficient at each sample size, serving as an indicator of network stability.  
* ar1_input_file: The input file for AR1 analysis, containing network topological parameters (e.g., node number, edge number, or density) calculated across the sample-size gradient.  
* ar1_windows: Defines the sliding window sizes (e.g.,15) for calculating rolling statistics like standard deviation or autocorrelation in the AR1 stability analysis.

### **2. Case Studies**
#### **2.1 Determination of N<sub>min</sub>**
In this part, the long-term time-series dataset from Station ALOHA is used as an example data set to determine the required N<sub>min</sub> for this dataset.
```r
library(SAFEMCN)

results <- analyze_network_topology(
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
  fit_formula = FALSE,
  fit_start_size = 11,
  predict_end_size = 58,
  run_ar1 = TRUE,
  ar1_input_file = "node_num_combined.csv",
  ar1_windows = c(15)
)
```

<p align="center">
  <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/ALOHA-test.png" 
       width="65%" 
       alt="Variations in network topological parameters and the lag-1 autocorrelation (AR1) coefficient" />
</p>


#### 📌 **Result Interpretation:**
The *SAFEMCN* package visualizes the variation of network topology parameters with sample size as a line graph, providing an intuitive representation of topological changes, while also outputting the corresponding topological parameter values for each sample size (e.g., node_num_original.csv, edge_num_original.csv) and saving them in a designated folder.   

Furthermore, the AR1 module of the *SAFEMCN*  package takes these output files (e.g., node_num_original.csv) as input to identify the Nmin threshold: the lowest AR1 value (marked by a rapid decrease followed by a rapid increase) is determined as Nmin, indicating the emergence of a topological plateau.   

Additionally, the AR1 module outputs a CSV file containing the AR1 coefficient for each sample size (e.g., node_num_original_window15_AR1.csv) and generates a scatter plot showing how the AR1 coefficient changes with sample size. Finally, the module also generates an AR1_min.csv file, which records the determined Nmin value and its corresponding AR1 value.

#### **2.2 Prediction of N<sub>min</sub>**
This part uses the ALOHA dataset as the sample data and employs the prediction module of *SAFEMCN* to predict Nmin.

<p align="center">
  <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/ALOHA-predicted.png" 
       width="65%" 
       alt="Variations in network topological parameters and the lag-1 autocorrelation (AR1) coefficient" />
</p>

#### 📌 **Result Interpretation:**
By extrapolating node and edge trajectories via the fitted exponential model (Figs 3a–b), the AR1-based N<sub>min</sub> was successfully predicted under limited sampling. Crucially, the predicted N<sub>min</sub> tightly matched the empirical thresholds derived from the complete 58-sample dataset (node: 25 vs. 27 samples; edge: 29 vs. 28 samples). This minimal deviation validates the framework's reliability. Furthermore, networks constructed with sample sizes exceeding these predicted thresholds exhibited convergent and stable topological properties, such as network density and average path length (Figs 3c–d).  

Consequently, combining mechanistic extrapolation with the AR1 stability criterion provides a robust, practical solution for optimized sampling design when large datasets are unfeasible.
